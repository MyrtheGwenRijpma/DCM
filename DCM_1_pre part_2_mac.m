%% Set directories
main_dir = '/shared/macdata/groups/rankin/rsfMRI_Library/';
ROIfolder = '/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/ROIs/';
glmpath = 'GLM_cos';

% Opening excel spreadsheet with PIDN, DCDate named 'Sample_NC_bvFTD_8mm_2020'
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the PIDN file');
MCINT = readtable(fullfile(PathName,Filename));

% formatting PIDN
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% formatting DCDATE
dcdate = MCINT.DCDate; %correct
date_input = datetime(dcdate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(date_input,formatOut);

% create subject variable
subj = strcat(pidn3,'_',DCDate);

%% Retrieve ROIs under the file name: 'ROIs_full_8mm'
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select your ROI file');
ROItable = readtable(fullfile(PathName,Filename));

for i=1:height(ROItable)
    rois{i,1} = [ROItable.ROI_name{i},'_',num2str(ROItable.size(i)),'mm'];
end

%----------------------------------------%

% Initialising SPM (add path under /Volumes/macdata/groups/rankin/spm12_new)
spm('Defaults','fMRI');
spm_jobman('initcfg');
 
%% Create VOIs via function style
spmthresh=0.05;

createVOI(subj, main_dir, glmpath, ROIfolder, rois, spmthresh);
