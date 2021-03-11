%% Set paths
addpath('/shared/macdata/groups/rankin/spm12_new')
addpath('/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/DCMscripts/functions')

%% Set directories and files
pathFolder = '/shared/macdata/groups/rankin/rsfMRI_Library';
dcm_model = 'SN_bilateral'; % all 10 SN ROIs (bilateral)
dcm_general = fullfile('/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/models/SN_bilateral/SN_level3');
dcm_results = fullfile(dcm_general,'DCM_results');

% subject DCM file
subj_dcm_file = 'SN_bilateral.mat'; %DCM for all 10 SN nodes selected (bilaterally)

%Check if DCM folder exist
if exist(dcm_folder,'dir') == 7 
    disp('Folder Exist!')
else % Create folder!
    mkdir(dcm_folder);
    disp('Folder created!')
end

%Check if DCM results folder exist
if exist(dcm_results,'dir') == 7 
    disp('Folder Exist!')
else % Create folder!
    mkdir(dcm_results);
    disp('Folder created!')
end

%----------------------------------------%

% Opening excel spreadsheet with PIDN, DCDate named 'Sample_NC_bvFTD_8mm_2020'
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the EXCEL file');
MCINT = readtable(fullfile(PathName,Filename));

% formatting PIDN
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% formatting DCDATE
date_input = datetime(MCINT.DCDate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(MCINT.DCDate,formatOut);

% create subject variable
subj_all = strcat(pidn3,'_',DCDate);
MCINT.PIDN_DCDate = subj_all;
clear pidn pidn2 pidn3 subj_all formatOut Filename PathName date_input DCDate

% filter and sort subjects
subjmatrix = [MCINT.PIDN_DCDate, MCINT.Dx];
[~,idx] = sort(subjmatrix(:,2));
nc_bvftd_matr_sort = subjmatrix(idx,:); % generating table of bvFTDs and NC subjects
subj = nc_bvftd_matr_sort(:,1);

% Retrieve the number of NCs and bvFTDs
n_nc = nnz(strcmp(nc_bvftd_matr_sort(:,2),'NC'));
n_bvftd = nnz(strcmp(nc_bvftd_matr_sort(:,2),'bvFTD'));

% Filter per DX group
nc_matr = nc_bvftd_matr_sort(strcmp(nc_bvftd_matr_sort(:,2),'NC'),:);
subj_NC = nc_matr(:,1);

bvFTD_matr = nc_bvftd_matr_sort(strcmp(nc_bvftd_matr_sort(:,2),'bvFTD'),:);
subj_bvFTD = bvFTD_matr(:,1);

%----------------------------------------%

%% Retrieve ROIs under the file name: 'ROIs_full_8mm'
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the ROI file');
ROItable = readtable(fullfile(PathName,Filename));

% Get list of all possible connections
counter = 1;
for i = 1:height(ROItable)
    for j = 1:height(ROItable)
        rois{1, counter+j-1} = [ROItable.ROI_name{i},'_',ROItable.ROI_name{j}];
    end
    counter = counter + 6;
end

%% Moving DCM_rest.mat files of interest into another folder
for ii=1:size(subj,1)
   
    % change this accordingly to _rest or _cort
    source = fullfile(pathFolder,subj{ii},'GLM_cos',subj_dcm_file); %change this accordingly
    destination = fullfile(dcm_folder,['DCM_',subj{ii},'_',subj_dcm_file]); %change this accordingly
    copyfile(source,destination);

end

%----------------------------------------%

%%first level PEB
% Collate DCMs into a GCM file
GCM_NC = fullfile(dcm_folder,strcat('DCM_',subj_NC,'_',subj_dcm_file));
GCM_bvFTD = fullfile(dcm_folder,strcat('DCM_',subj_bvFTD,'_',subj_dcm_file));


%% Estimate a first level PEB
tic; [GCM_NC, DCM_NC] = spm_dcm_peb_fit(GCM_NC); toc;
tic; [GCM_bvFTD, DCM_bvFTD] = spm_dcm_peb_fit(GCM_bvFTD); toc;


% Move DCM results to DCM results folder
save(fullfile(dcm_results,'GCM_DCM_fit_NC.mat'), 'GCM_NC');
save(fullfile(dcm_results,'GCM_DCM_fit_bvFTD.mat'), 'GCM_bvFTD');

%% Put first level DCM parameters into an excel spreadsheet (correlation and prediction)
[subCM_NC, subPM_NC] = getCMPM(GCM_NC, nc_matr, rois);
[subCM_bvFTD, subPM_bvFTD] = getCMPM(GCM_bvFTD, bvFTD_matr, rois);
%[subCM_NC, subPM_NC] = getCMPM(GCM_NC, nc_matr(:,1:2), rois);

% Save to excel for NC
writetable(subCM_NC,fullfile(dcm_results,'SN_bi_subCM_NC.csv'));
writetable(subPM_NC,fullfile(dcm_results,'SN_bi_subPM_NC.csv'));

% Save to excel for bvFTD
writetable(subCM_bvFTD,fullfile(dcm_results,'SN_bi_subCM_NC.csv'));
writetable(subPM_bvFTD,fullfile(dcm_results,'SN_bi_subPM_NC.csv'));

%----------------------------------------%

%% Second level PEB
% Specify PEB model settings
M = []; 
M = struct();
M.alpha = 1;
M.beta = 16;
M.hE = 0;
M.hc = 1/16;
M.Q = 'single';

% Specifiy the design matrix. It should start with a constant column (1).
NC = str2double(nc_matr(:,3));
bvFTD = str2double(bvFTD_matr(:,3)); 

% grand mean NC
M.X = [ones(n_nc,1), ...  % grand mean
    
% Choose a field
field = {'A'};

% Estimate second level PEB for NC
PEB_NC = spm_dcm_peb(GCM_NC,M,field);
save(fullfile(dcm_results,'PEB2_NC.mat'), 'PEB_NC');

% grand mean bvFTD
M.X = [ones(n_bvftd,1), ...  % grand mean

% Estimate second level PEB for bvFTD
PEB_bvFTD = spm_dcm_peb(GCM_bvFTD,M,field);
save(fullfile(dcm_results,'PEB2_bvFTD.mat'), 'PEB_bvFTD');

% Bayesian Model Averaging (both DX groups)
BMA_NC = spm_dcm_peb_bmc(PEB_NC);
BMA_bvFTD = spm_dcm_peb_bmc(PEB_bvFTD);
save(fullfile(dcm_results,'BMA_NC.mat'), 'BMA_NC');
save(fullfile(dcm_results,'BMA_bvFTD.mat'), 'BMA_bvFTD');

% Review results for PEB and BMA
spm_dcm_peb_review(BMA_NC,GCM_NC);
spm_dcm_peb_review(BMA_bvFTD,GCM_bvFTD);

% Perform leave-one-out cross validation (GCM,M,field are as before)
class=spm_dcm_loo(GCM,M,field);

%----------------------------------------%

%% Third Level PEB
% specificy design matrix for NC and bvFTD combined (third PEB)
X3 = [ones(n_nc+n_bvftd,1),cat(1,ones(n_nc,1)*-1,ones(n_bvftd,1)*1)];
M.X = X3;
PEBs = {PEB_NC; PEB_bvFTD};
PEB3 = spm_dcm_peb(PEBs,X3);

%Bayesian Model Averaging for third PEB
BMA_difference = spm_dcm_peb(PEB3);

%GCM file for NCa nd bvFTD combined
GCMs = [GCM_NC; GCM_bvFTD];
spm_dcm_peb_review(PEB3, GCMs);

% save PEB3 variables;
save(fullfile(dcm_results,'GCM_DCM_fit_PEB3.mat'), 'GCMs');
save(fullfile(dcm_results,'PEB_DCM_fit_PEB3.mat'), 'PEBs');

%% Get group parameters
totalparams = length(rois)*2;
BMA_Ep = BMA.Ep(length(rois)+1:totalparams);
BMA_Ep = full(BMA_Ep);
BMA_Pp = BMA.Pp(length(rois)+1:totalparams);
BMA_Pp = full(BMA_Pp);

PEB_Ep = PEB.Ep(length(rois)+1:totalparams);
PEB_Ep = full(PEB_Ep);
PEB_Pp = PEB.Pp(length(rois)+1:totalparams);
PEB_Pp = full(PEB_Pp);

CM = num2cell(CM);
PM = num2cell(PM);
