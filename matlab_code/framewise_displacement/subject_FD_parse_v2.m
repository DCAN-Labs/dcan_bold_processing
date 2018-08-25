function subject_FD_parse_v2(FD_movement_files,skip_seconds,epi_TR,brain_radius_in_mm,result_dir)
% Code from OHSU Fair Neuroimaging Lab
% Written (2014) and commented (9/8/2015) by Eric Earl
% 
% FD parsing Nipype wrapper pulling in list of motion files
% 
% Modified to read (and ignore) DVAR_pre_reg, DVAR_post_reg, and DVAR_post_all
% All variables (skips, TR, CurrentDirectory) regularly set within Nipype
% 
% Iterates over range of values from 0.1 to maximum of 0.5 FD threshold
% Outputs to .mat file
% 
% Examples:
% FD_movement_files = {'motion_numbers0.txt', 'motion_numbers1.txt', 'motion_numbers2.txt'};
% FD_movement_files can be FNL TXT files or 4dfp DAT files
% varargin (brain_radius_in_mm) only needs to be provided for DAT file inputs
% 
% NOTE: need to comment on how FD_movement_file is structured

% Calculates the common range of FD's from 0 to 0.5 in steps of 0.01
FD_range = 0:0.01:0.5; 
FD_data = {};
run_count = length(FD_movement_files); % each movement file is one run

for i = 1:length(FD_range) % for each FD threshold
    
    FD_data{i}.skip = ceil(skip_seconds/epi_TR); % Frames to skip
    FD_data{i}.epi_TR = epi_TR;   % EPI/BOLD/Functional Run Repetition Time
    FD_data{i}.FD_threshold = FD_range(i); % FD Threshold
    FD_data{i}.frame_removal = []; % Initialize empty frame_removal vector
    
    FD_all_runs = [];
    
    for j = 1:run_count % for each run/series
        movement_file  = FD_movement_files{j};
        
        % Pad FD vector with a false frame (difference of frame "0" and 1) 
        % before the first actual FD number (the difference of frame 1
        % and 2)
        % ADD ANY OTHER PARSER FOR ANY OTHER FILE FORMAT HERE
        switch movement_file(end-3:end)
            case '.txt' % Read FD values straight from FNL motion numbers files
                fid = fopen(movement_file);
                FD_cell = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %*s %*s','HeaderLines',1);% Read FD values from motion numbers files
                fclose(fid);
                FD_vector = cat(1,0,FD_cell{1});
            case '.dat' % Calculate FD valuescondor from 4dfp DAT files
                FD_vector = cat(1,0,calc_FD(movement_file,brain_radius_in_mm));
            otherwise
                disp('Input is neither a .dat or .txt file.  Exiting.')
                return
        end
        
        % The removal vector is calculated by concatenating the removal
        % vector from the previous run/series iteration to the current
        % run/series iteration's removal vector coming out of the
        % pre_post_motion_censoring.m function
        FD_data{i}.frame_removal = cat(1,FD_data{i}.frame_removal, ...
            pre_post_motion_censoring(FD_vector,FD_data{i}.FD_threshold,FD_data{i}.skip));
        
        % Make an FD for all runs vector as well
        FD_all_runs = cat(1,FD_all_runs,FD_vector);
        
    end
    
    FD_data{i}.format_string = format_generator(FD_data{i}.frame_removal); % output a fromat string for 4dfp usage
    FD_data{i}.total_frame_count = length(FD_data{i}.frame_removal); % output the total frame count
    FD_data{i}.remaining_frame_count = FD_data{i}.total_frame_count - sum(FD_data{i}.frame_removal); % output the remaining frame count
    FD_data{i}.remaining_seconds = FD_data{i}.remaining_frame_count * FD_data{i}.epi_TR; % convert the remaining frame count to seconds using TR
    FD_data{i}.remaining_frame_mean_FD = mean(FD_all_runs(FD_data{i}.frame_removal==0)); % calculate the mean FD of remaining frames
    
end

% Save output FD_data cell array of structs to an FD.mat file in the
% present working directory
% this_path=pwd;
%save([this_path filesep 'FD.mat'],'FD_data');
save([result_dir filesep 'FD.mat'],'FD_data');
