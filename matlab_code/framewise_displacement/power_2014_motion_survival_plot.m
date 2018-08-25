function [subject_list,subject_survival] = power_2014_motion_survival_plot(group_folder_path,subject_list_file,FD_threshold)
%Displays a colour map of the number of seconds remaining at various FD thresholds

% Inputs:
% group_folder_path     Full filepath to where your subject visit folders
%                       are located with timecourses and FD mat files.
% subject_list_file     The full path to your subject list .txt file
% 
% Outputs:
% subject_list          A cell of subject names
% 

% read in subject list
fid=fopen(subject_list_file);
vc_names=textscan(fid,'%s');
fclose(fid);

subject_list = vc_names{1};
subject_count = length(subject_list);

subject_survival = [];

% loop through subjects
for i = 1:subject_count % for each subject
    sub_name = subject_list{i};
    disp(['Starting analyzing visit ' sub_name])
    subject_dir = [group_folder_path filesep sub_name];
    
    FD_path = [subject_dir filesep 'motion' filesep 'FD.mat'];
    
    try
        % for every subject get remaining seconds at the FD threshold
        subject_survival = cat(1,subject_survival,read_FD(FD_path,FD_threshold,'remaining_seconds')/60);
    catch FD_mat_exception
        disp(['Excluding ' sub_name ' (Missing or broken FD.mat link): ' FD_path])
    end
    
end

