%% Set directories for datasets
main_dir = '/shared/macdata/groups/rankin/rsfMRI_Library/';
glmpath = 'GLM_cos';

% Opening excel spreadsheet with PIDN, DCDate named 'Sample_NC_bvFTD_8mm_2020'
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the PIDN DCDate EXCEL file');
MCINT = readtable(fullfile(PathName,Filename));

% formatting PIDN
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% formatting DCDATE
date_input = datetime(MCINT.DCDate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(date_input,formatOut);

% create subject variable
subj = strcat(pidn3,'_',DCDate);

%----------------------------------------%

% Set parameters
TR = 2; 
n_vol = 235;

% Initialising SPM (add path under /Volumes/macdata/groups/rankin/spm12_new)
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Create DCM regressor
DCMregressor(subj, main_dir, glmpath, TR, n_vol);
