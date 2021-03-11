%% Set directories
main_dir='/shared/macdata/groups/rankin/rsfMRI_Library/';

%% Set up filenames
DCMfilename = 'SN_bilateral.mat'; %DCM for all 10 SN nodes selected (bilaterally)
DCM_estimateoutput = ['estim_output_',DCMfilename];
DCM_estimate = ['estimate_',DCMfilename];
DCM_filenames = {DCMfilename, DCM_estimateoutput, DCM_estimate};

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
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the ROI file');
ROItable = readtable(fullfile(PathName,Filename));

for i=1:height(ROItable)
    rois{i,1} = ['VOI_',ROItable.ROI_name{i},'_',num2str(ROItable.size(i)),'mm_1.mat'];
end

%----------------------------------------%

% Initialising SPM (add path under /Volumes/macdata/groups/rankin/spm12_new)
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Specifying the DCM
sTable = specifyDCM(subj, main_dir, DCM_filenames, rois)
