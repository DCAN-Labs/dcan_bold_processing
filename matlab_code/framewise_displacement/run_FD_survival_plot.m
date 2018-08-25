% Script to take FD_survival_plot input directory and path and input them
% to the function as well as modify a few figure fields post-function run

%% Example argument formatting:

% group_folder_path = '/group_shares/FAIR_LAB2/CYA/analyses_v2/Music-HumanAdult-Auckland_OHSU/GEN3-SEEDCORREL_1';
% subject_list_file = '/group_shares/FAIR_LAB2/CYA/analyses_v2/Music-HumanAdult-Auckland_OHSU/GEN3-SEEDCORREL_1/subjectlist.txt';
%% 7/16/2014 LKV added :
% output of the subject_survival matrix to the workspace
% output of a list of subjects whose FD data is missing
%% General solution asking user to navigate to group folder and subject list

group_folder_path = uigetdir('/group_shares/FAIR_LAB2/CYA/analyses_v2','Select group folder path');
[subjlistfile,subjlistpath] = uigetfile([group_folder_path filesep '*.txt'],'Select subject list file');
subject_list_file = [subjlistpath filesep subjlistfile];

%% Call function using user inputs from this script

[subject_list_cell, subject_survival, missing_data] = FD_survival_plot(group_folder_path,subject_list_file);

%% Changing figure properties

set(gca,'ytick',1:length(subject_list_cell),'yticklabel',subject_list_cell);
ylabel('Subject Visits')
title('Remaining time in minutes after FD thresholding')
