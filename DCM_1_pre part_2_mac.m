%% Set directories
main_dir = '/shared/macdata/groups/rankin/rsfMRI_Library/';
ROIfolder = '/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/ROIs/';
glmpath = 'GLM_cos';

%% Get Subjects list
% Opening File directory
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the PIDN file');
MCINT = readtable(fullfile(PathName,Filename));

% Create PIDN_DCDATE
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
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select your ROI file');
ROItable = readtable(fullfile(PathName,Filename));

for i=1:height(ROItable)
    rois{i,1} = [ROItable.ROI_name{i},'_',num2str(ROItable.size(i)),'mm'];
end


%% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
 
%% Create VOIs via function style
spmthresh=0.05;

createVOI(subj, main_dir, glmpath, ROIfolder, rois, spmthresh);



%% Looping Creation of VOIs
% for ii=1:size(subj,1)
% 
%     dir_ind = fullfile(main_dir, subj{ii});
%     GLM_dir = fullfile(dir_ind,'GLM');
%     
%     %Create VOI folder to store the VOIs
%     mkdir(GLM_dir,'VOI');
%     
%     VOI_dir = fullfile(GLM_dir,'VOI');
%     
%     %---------------------%
%     % VOLUMES OF INTEREST %
%     %---------------------%
%     
%     clear matlabbatch SPM  
%     
%     for jj=1:size(rois)
%         
%         %--------------------------------------------------%
%         % EXTRACTING TIME SERIES FOR ALL ROIS IN EACH PIDN %
%         %--------------------------------------------------%
%         
%         matlabbatch{jj}.spm.util.voi.spmmat = cellstr(fullfile(dir_ind,'GLM','SPM.mat'));
%         matlabbatch{jj}.spm.util.voi.adjust = 0; 
%         matlabbatch{jj}.spm.util.voi.session = 1; % session 1
%         matlabbatch{jj}.spm.util.voi.name = rois{jj};
%         matlabbatch{jj}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
%         matlabbatch{jj}.spm.util.voi.roi{1}.spm.contrast = 1;  % F test
%         matlabbatch{jj}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
%         matlabbatch{jj}.spm.util.voi.roi{1}.spm.thresh = 0.05;
%         matlabbatch{jj}.spm.util.voi.roi{1}.spm.extent = 0;
%         %matlabbatch{jj}.spm.util.voi.roi{2}.spm.mask.contrast = 1; 
%         matlabbatch{jj}.spm.util.voi.roi{2}.mask.image = cellstr(fullfile(ROIfolder,[rois{jj},'.nii,1']));
% 
%         matlabbatch{jj}.spm.util.voi.roi{2}.mask.threshold = 0.5;
%         matlabbatch{jj}.spm.util.voi.expression = 'i1 & i2';
%     end
%     
%     save(fullfile(VOI_dir,'VOI_extraction.mat'),'matlabbatch');
%     spm_jobman('run',matlabbatch);
%     
% end

