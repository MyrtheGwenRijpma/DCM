%% Add paths SPM12_new and funtions
addpath('/smb://MAC-SRV-NAS2.UCSF.EDU/macdata/groups/rankin/spm12_new')
addpath('/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/DCMscripts/functions')

%% Set directories and files
% Change the folder names according to the model
pathFolder = '/shared/macdata/groups/rankin/rsfMRI_Library/';
dcm_model = 'SN_bilateral'; % This is where we need to change name according to specified model
dcm_general = fullfile('/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/models/SN_bilateral/SN_level3');
dcm_folder = fullfile(dcm_general,'DCM_files');
dcm_results = fullfile(dcm_general,'DCM_results'); %change according to covariates and subject#
PEB3_results = fullfile('Users/mrijpma/Desktop/barchart_DCM');

% subject DCM file
subj_dcm_file = 'SN_bilateral.mat';

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

% List all directories (correspond to the subject names)

%% Get PIDN_DCDate

% Opening File directory
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the EXCEL file');
MCINT = readtable(fullfile(PathName,Filename));

% Creating PIDN_DCDATE
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% changing format of the date
date_input = datetime(MCINT.DCDate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(MCINT.DCDate,formatOut);

%PIDNs = MCINT.PIDN_DCDATE; % PIDN_DCDATE variable
subj_all = strcat(pidn3,'_',DCDate);

MCINT.PIDN_DCDate = subj_all;
% clear unwanted variables
clear pidn pidn2 pidn3 subj_all formatOut Filename PathName date_input DCDate

%% SUBJECT FILTER, AND SORT
subjmatrix = [MCINT.PIDN_DCDate, MCINT.Dx, cellstr(num2str(MCINT.RSMS))]; % RSMS as covariate

% Sort
[~,idx] = sort(subjmatrix(:,2));
nc_bvftd_matr_sort = subjmatrix(idx,:); % Gets table

subj = nc_bvftd_matr_sort(:,1);

% Get number of NCs and bvFTDs
n_nc = nnz(strcmp(nc_bvftd_matr_sort(:,2),'NC'));
n_bvftd = nnz(strcmp(nc_bvftd_matr_sort(:,2),'bvFTD'));

% Filter for each Dx group
nc_matr = nc_bvftd_matr_sort(strcmp(nc_bvftd_matr_sort(:,2),'NC'),:);
subj_NC = nc_matr(:,1);

bvFTD_matr = nc_bvftd_matr_sort(strcmp(nc_bvftd_matr_sort(:,2),'bvFTD'),:);
subj_bvFTD = bvFTD_matr(:,1);


%% Get ROIs list
% Opening File directory
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

%% Moving interested DCM_rest.mat into another folder

% I copy all the DCMs I am interested in a different folder

for ii=1:size(subj,1)
   
    % change this accordingly to _rest or _cort
    source = fullfile(pathFolder,subj{ii},'GLM_cos',subj_dcm_file); %change this accordingly
    destination = fullfile(dcm_folder,['DCM_',subj{ii},'_',subj_dcm_file]); %change this accordingly
    copyfile(source,destination);

end

% Collate DCMs into a GCM file
GCM_NC = fullfile(dcm_folder,strcat('DCM_',subj_NC,'_',subj_dcm_file));
GCM_bvFTD = fullfile(dcm_folder,strcat('DCM_',subj_bvFTD,'_',subj_dcm_file));

%% Estimate a first level PEB

% Fully estimate model
% for ii=1:size(GCM,1)
%     GCM(:,ii) = spm_dcm_fit(GCM(:,ii)); 
% end

% Fully estimate model
% use_parfor = true ;
% GCM = spm_dcm_fit ( GCM , use_parfor ) ;


% Use Bayesian Model Reduction to rapidly estimated DCMs 2-N for each subject if applicable
% if size(GCM,2) > 1
%    GCM = spm_dcm_bmr(GCM);
% end

% alternate between estimating DCMs and estimating group effects. 
% This is slower, but can draw subjects out of local optima towards the group mean.

tic; [GCM_NC, DCM_NC] = spm_dcm_peb_fit(GCM_NC); toc;
tic; [GCM_bvFTD, DCM_bvFTD] = spm_dcm_peb_fit(GCM_bvFTD); toc;


% Move DCM results to DCM results folder
save(fullfile(dcm_results,'GCM_DCM_fit_NC.mat'), 'GCM_NC');
save(fullfile(dcm_results,'GCM_DCM_fit_bvFTD.mat'), 'GCM_bvFTD');

%% Put first level DCM parameters into an excel spreadsheet
% for i = 1:size(subj,1)
%     subCM(i,:) = GCM{i}.Ep.A(:)'; % Connectivity matrix
%     subPM(i,:) = GCM{i}.Pp.A(:)'; % Probability matrix
% end
% 
% subCM = num2cell(subCM);
% subPM = num2cell(subPM);
% 
% % Combine matrix and demographics
% subCM_table = [nc_bvftd_matr_sort, subCM];
% subPM_table = [nc_bvftd_matr_sort, subPM];
% subCM_table = cell2table(subCM_table, 'VariableNames',['PIDN_DCDate' 'Dx' rois]);
% subPM_table = cell2table(subPM_table, 'VariableNames',['PIDN_DCDate' 'Dx' rois]);

[subCM_NC, subPM_NC] = getCMPM(GCM_NC, nc_matr, rois);
[subCM_bvFTD, subPM_bvFTD] = getCMPM(GCM_bvFTD, bvFTD_matr, rois);

[subCM_NC, subPM_NC] = getCMPM(GCM_NC, nc_matr(:,1:2), rois);

% Save to excel for NC
writetable(subCM_NC,fullfile(dcm_results,'SN_right_hemi_subCM_NC.csv'));
writetable(subPM_NC,fullfile(dcm_results,'SN_right_hemi_subPM_NC.csv'));

% Save for bvFTD
writetable(subCM_bvFTD,fullfile(dcm_results,'SN_right_hemi_subCM_NC.csv'));
writetable(subPM_bvFTD,fullfile(dcm_results,'SN_right_hemi_subPM_NC.csv'));


[subCM, subPM] = getCMPM(GCMs, subjmatrix(:,1:2), rois);
writetable(subCM,fullfile(dcm_results,'SN_right_hemi_subCM3.csv'));
writetable(subPM,fullfile(dcm_results,'SN_right_hemi_subPM3.csv'));



%% Estimate a second level PEB (Parametric Empirical Bayes) model - include RSMS
% Specify PEB model settings
M = []; 
M = struct();
M.alpha = 1;
M.beta = 16;
M.hE = 0;
M.hc = 1/16;
M.Q = 'single';

% Specifiy the design matrix. It should start with a constant column (1).

% Create design matrix
RSMS_NC = str2double(nc_matr(:,3)); %define RSMS: e.g. [67;50;51]
RSMS_bvFTD = str2double(bvFTD_matr(:,3)); %define RSMS: e.g. [67;50;51]

% RSMS main regressor for NC
M.X = [ones(n_nc,1), ...  % grand mean
    RSMS_NC];      % RSMS regressor
    
% Choose a field
field = {'A'};

% Estimate second level PEB for NC
PEB_NC = spm_dcm_peb(GCM_NC,M,field);
save(fullfile(dcm_results,'PEB_RSMS_NC.mat'), 'PEB_NC');

% RSMS main regressor for bvFTD
M.X = [ones(n_bvftd,1), ...  % grand mean
    RSMS_bvFTD];      % RSMS regressor


% Estimate second level PEB for bvFTD
PEB_bvFTD = spm_dcm_peb(GCM_bvFTD,M,field);
save(fullfile(dcm_results,'PEB_RSMS_bvFTD.mat'), 'PEB_bvFTD');


% Compare nested PEB models. Decide which connections to switch off based on the 
% structure of each DCM for subject 1.
% BMA = spm_dcm_peb_bmc(PEB(1), GCM(1,:));

% Search over nested models
% Rather than compare specific hypotheses, 
% you may wish to simply prune away any parameters from the PEB which don't contribute to the model evidence.
BMA_NC = spm_dcm_peb_bmc(PEB_NC);
BMA_bvFTD = spm_dcm_peb_bmc(PEB_bvFTD);
save(fullfile(dcm_results,'BMA_RSMS_NC.mat'), 'BMA_NC');
save(fullfile(dcm_results,'BMA_RSMS_bvFTD.mat'), 'BMA_bvFTD');

% Review results for PEB and BMA
spm_dcm_peb_review(BMA_NC,GCM_NC);
spm_dcm_peb_review(BMA_bvFTD,GCM_bvFTD);

% Perform leave-one-out cross validation (GCM,M,field are as before)
class=spm_dcm_loo(GCM,M,field);

%% Third Level PEB

X3 = [ones(n_nc+n_bvftd,1),cat(1,ones(n_nc,1)*-1,ones(n_bvftd,1)*1)];
M.X = X3;
PEBs = {PEB_NC; PEB_bvFTD};
PEB3 = spm_dcm_peb(PEBs,X3);

% BMA_difference = spm_dcm_peb(PEB3);

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