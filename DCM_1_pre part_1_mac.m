%% Set directories
% cd '/shared/macdata/groups/rankin/Users/Gianina/rsFMRI/DCM/'; % change directory
main_dir = '/shared/macdata/groups/rankin/rsfMRI_Library/';
glmpath = 'GLM_cos';

%% Get Subjects list
% Opening File directory 
% Excel spreadsheet with PIDN, DCDate
[Filename,PathName] = uigetfile({'*.xlsx';'*.xls';'*.*'},'Select the EXCEL file');
MCINT = readtable(fullfile(PathName,Filename));

% Create PIDN_DCDATE
pidn = num2str(MCINT.PIDN);
pidn2 = cellstr(pidn);
pidn3 = regexprep(pidn2, '\W', '');

% changing format of the date
date_input = datetime(MCINT.DCDate,'ConvertFrom','excel');
formatOut = 'mmddyyyy';
DCDate = datestr(date_input,formatOut);

%PIDNs = MCINT.PIDN_DCDATE; % PIDN_DCDATE variable
subj = strcat(pidn3,'_',DCDate);

%% Set parameters
TR = 2; 
n_vol = 235; % specify accordingly 

%% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');


%% Create DCM regressor
DCMregressor(subj, main_dir, glmpath, TR, n_vol);

%% Start looping
% for ii=1:size(subj,1)
%     
%     % Print subject
%     disp(['creating regressor file for ', subj{ii}]);
%     
%     %--------------------------------%
%     % create multiple regressor file %
%     %--------------------------------%
% 
%     % Retrieve motion and regressor files
%     motionfile = 'rsfmri/interfmri_TRCNnSFmDI/motion_corr/rp_avol_006.txt';
%     regressorfile = '/rsfmri/processedfmri_TRCNnSFmDI/stats_FC_R_vAI_SN_sphere_roi/timeseries/seed_nuisance_regressors.txt';
% 
%     % Load motion and regressor files
%     load(fullfile(main_dir,subj{ii},motionfile));
%     load(fullfile(main_dir,subj{ii},regressorfile));
%     covar=[rp_avol_006,seed_nuisance_regressors(:,2:3)];
% 
%     % Save covariates
%     % eval(['save ',main_dir,subj{ii},'/covar.mat',' covar']);
% 
%     % DCT.  Creates discrete cosine set with frequencies ranging from the UL to the
%     % LL (default UL = 0.1Hz, LL = 1/128hz). Inputs are:
%     % dir = path for output SPM
%     % n = number of scans
%     % covar =  multiple regressor file
%     %--------------------------------------------------------------------------
% 
%     dir_ind = fullfile(main_dir, subj{ii});
%     
%     % specify accordingly
%     TR=2; 
%     n_vol=175; 
%     
%     [n_cols,R] = spm_glm_rest_dct(dir_ind,TR,n_vol,covar); % no mat file given. output is a text file
%     
%     % To Do: Update directory
%     save(fullfile(dir_ind,'glm_regr.mat'), 'R');
% 
%     % -----------------------------------------------------
%     % Extra steps
%     % load ([main_dir,subj{ii},'\glm_regr.mat']); 
%     % mreg_full_dct=[mreg_full_dct];
%     % save (strcat(main_dir,subj{ii},'\mreg_full_dct.txt'),'mreg_full_dct')
%     % -----------------------------------------------------
% 
%     %eval(['!gunzip ', dir_ind, '/',subj{ii},'_rfMRI_2_TT_N27_NonLinear_Warping_smoothed_GMS_masked.nii.gz'])
% 
%     %-------------------------------------------%
%     % GLM SPECIFICATION, ESTIMATION & INFERENCE %
%     %-------------------------------------------%
% 
%     % factors = load(fullfile(dir_ind,'glm_regr.mat')); % NEED TO CHECK!!!
%     factors = load(fullfile(dir_ind,'glm_regr.mat')); % adapted
%     f = spm_select('FPList', fullfile(dir_ind,'rsfmri/interfmri_TRCNnSFmDI/images/'), 'swuavol_.*');
%     % f = cellstr(f)
%     
%     % img = cellstr(strcat(f, ',1'));
%     for jj=1:n_vol
%     img{jj,1}=[f(jj,:),',',int2str(1)];
%     end
% 
%     clear matlabbatch
% 
%     % OUTPUT DIRECTORY
%     %--------------------------------------------------------------------------
%     matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(dir_ind);
%     matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';
% 
% 
%     % MODEL SPECIFICATION
%     %--------------------------------------------------------------------------
%     matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(dir_ind,'GLM'));
%     matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
%     matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 2;
%     matlabbatch{2}.spm.stats.fmri_spec.sess.scans = img(:,1);
%     matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 100;
%     matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = cellstr([main_dir, subj{ii},'/glm_regr.mat']);
%     matlabbatch{1,2}.spm.stats.fmri_spec.bases  = struct('none',1);
% 
%     % MODEL ESTIMATION
%     %--------------------------------------------------------------------------
%     matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(dir_ind,'GLM','SPM.mat'));
% 
%     % INFERENCE
%     %--------------------------------------------------------------------------
%     matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(dir_ind,'GLM','SPM.mat'));
%     matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
%     matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(n_cols);
% 
%     save(fullfile(subj_DCM_fold,'model_spec_inference.mat'), 'matlabbatch');
%     spm_jobman('run',matlabbatch);
%     
%     disp(['creating regressor file for ', subj{ii}, ' complete!!']);
%     
%     clear matlabbatch
%     clear matlabbatch SPM
% 
% end
