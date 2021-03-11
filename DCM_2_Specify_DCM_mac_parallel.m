%cd 'R:\groups\rankin\rsfMRI_Library\';

pathFolder='/shared/macdata/groups/rankin/rsfMRI_Library/';
main_dir='/shared/macdata/groups/rankin/rsfMRI_Library/';
%pathFolder='D:\Winson\UCSF\projects\DCM\scans\';

%% Set up filenames
DCMfilename = 'SN_bilateral.mat'; %change name depending on model
DCM_estimateoutput = ['estim_output_',DCMfilename];
DCM_estimate = ['estimate_',DCMfilename];
DCM_filenames = {DCMfilename, DCM_estimateoutput, DCM_estimate};

%% Get Subjects list

% Opening File directory
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the PIDN file');
MCINT = readtable(fullfile(PathName,Filename));

% Creating PIDN_DCDATE
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% changing format of the date
dcdate = MCINT.DCDate; %correct
date_input = datetime(dcdate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(date_input,formatOut);

%PIDNs = MCINT.PIDN_DCDATE; % PIDN_DCDATE variable
subj = strcat(pidn3,'_',DCDate);


%% Retrieve ROIs
% Opening File directory
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the ROI file');
ROItable = readtable(fullfile(PathName,Filename));

for i=1:height(ROItable)
    rois{i,1} = ['VOI_',ROItable.ROI_name{i},'_',num2str(ROItable.size(i)),'mm_1.mat'];
end


%% SPM

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

sTable = specifyDCM(subj, main_dir, DCM_filenames, rois)
