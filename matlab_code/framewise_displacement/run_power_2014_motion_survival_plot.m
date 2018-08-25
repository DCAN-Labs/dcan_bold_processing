% Script to take FD_survival_plot input directory and path and input them
% to the function as well as modify a few figure fields post-function run
% 

%% Example argument formatting:

% group_folder_path = '/group_shares/FAIR_LAB2/CYA/analyses_v2/Music-HumanAdult-Auckland_OHSU/GEN3-SEEDCORREL_1';
% subject_list_file = '/group_shares/FAIR_LAB2/CYA/analyses_v2/Music-HumanAdult-Auckland_OHSU/GEN3-SEEDCORREL_1/subjectlist.txt';

%% General solution asking user to navigate to group folder and subject list

group_folder_path = uigetdir('/group_shares/FAIR_LAB2/CYA/analyses_v2','Select group folder path');
[subjlistfile,subjlistpath] = uigetfile([group_folder_path filesep '*.txt'],'Select subject list file');
subject_list_file = [subjlistpath filesep subjlistfile];
dlgcell = inputdlg({'FD Threshold (0 to 0.5)'},'Enter FD Threshold',1,{'0.2'});
FD_threshold = str2double(dlgcell{1});

%% Call function using user inputs from this script

subject_list_cell = power_2014_motion_survival_plot(group_folder_path,subject_list_file,FD_threshold);

%% Changing figure properties

set(gca,'ytick',1:length(subject_list_cell),'yticklabel',subject_list_cell);
ylabel('Subject Visits')
title('Remaining time in minutes after FD thresholding')
