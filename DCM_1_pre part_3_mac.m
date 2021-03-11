%% Set directories
main_dir = '/shared/macdata/groups/rankin/rsfMRI_Library/';
ROIfolder = '/shared/macdata/groups/rankin/Users/Myrthe/2019/rSMS_DCM/RSMS-DCM/ROIs/';

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
    rois{i,1} = ['VOI_',ROItable.ROI_name{i},'_',num2str(ROItable.size(i)),'mm_1.mat']
end

%% Move VOIs to VOI folder

% Get VOIs to move
VOIs = '*mm*'; 
transferVOI(subj, main_dir, VOIs)
% for iii=1:size(subj)
%     
%     % Set folder
%     ind_dir = fullfile(main_dir,subj{iii},'GLM_cos');
%     VOIfold = fullfile(ind_dir,'VOI');
%     
%     % Move files to folder
%     try
%         movefile(fullfile(ind_dir, VOIs),fullfile(VOIfold));
%         msg2 = [num2str(iii), '. Files for ', subj{iii},' moved to folder: ',voifolder];
%         disp(msg2);
%     catch
%         msg2 = [num2str(iii), '. No files found for ', subj{iii}];
%         disp(msg2);  
%     end
% end

%% Loop to create the time series figures
% for ii=1:size(subj)
%     
%     ind_dir = fullfile(main_dir,subj{ii},'GLM_cos', 'VOI');
%     
%     %% Load the VOIs
%     for jj=1:size(rois)
%         try
%             VOImat = load (fullfile(ind_dir,rois{jj}));
%             Y = VOImat.Y;
%             xY = VOImat.xY;
% 
%             %Rename VOI to look neat
%             VOItitle = strrep(rois{jj},'_',' ');
%             VOItitle = VOItitle(1:end-6);
% 
%             %% Plot and save the Time series
%             h(1) = figure;
%             plot(Y);
%             axis([1 235 min(Y) max(Y)]);
%             xlabel('Scan number'), ylabel('Frequency'); 
%             title(['Time Series of ',VOItitle]);
% 
%             % Save figure as png file
%             nametosave = fullfile(ind_dir,[VOItitle,'_TS.png']);
%             saveas(h,nametosave);
% 
%             % Close the figure (NEED SAVE BEFORE CLOSING)
%             close(h);
%         catch
%             disp('Error! No VOI found')
%         end % End Try
%     end % End VOIs loop  
% end % End Subject loop