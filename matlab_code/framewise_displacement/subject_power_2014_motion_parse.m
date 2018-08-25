function subject_power_2014_motion_parse(movement_files,skip_seconds,epi_TR,expected_contiguous_frame_count)
% FD parsing Nipype wrapper pulling in list of motion files
% All variables (skips, TR, CurrentDirectory) set within Nipype
% Iterates over range of values from 0.1 to maximum of 0.5 FD threshold
% Outputs to .mat file
% FD_movement_files = {'motion_numbers0.txt', 'motion_numbers1.txt', 'motion_numbers2.txt'};
% FD_movement_files can be TXT or DAT files
% varargin (brain_radius_in_mm) only needs to be provided for DAT file inputs
% NOTE: need to comment on how FD_movement_file is structured

FD_range = 0:0.01:0.5;
DVAR_range = 0:1:50;
motion_data = cell(length(FD_range),length(DVAR_range));
run_count = size(movement_files,2);

for i = 1:length(FD_range)
    for k = 1:length(DVAR_range)
        
        motion_data{i,k}.skip = ceil(skip_seconds/epi_TR); % Frames to skip
        motion_data{i,k}.epi_TR = epi_TR;   % EPI/BOLD/Functional Run Repetition Time
        motion_data{i,k}.FD_threshold = FD_range(i); % FD Threshold
        motion_data{i,k}.DVAR_threshold = DVAR_range(k); % DVAR Threshold
        motion_data{i,k}.frame_removal = []; % Initialize empty frame_removal vector
        
        FD_all_runs = [];
        DVAR_all_runs = [];
        
        for j = 1:run_count
            movement_file  = movement_files{j};
            
            % Pad FD vector with a false frame (difference of frame "0" and 1)
            % before the first actual FD number (the difference of frame 1
            % and 2)
            switch movement_file(end-3:end)
                case '.txt' % Read FD values from motion numbers files
                    fid = fopen(movement_file);
                    motion_cell = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f','HeaderLines',1);% Read FD values from motion numbers files
                    fclose(fid);
                    FD_vector = cat(1,0,motion_cell{1});
                    DVAR_vector = cat(1,0,motion_cell{2});
                otherwise
                    disp('Input is not a .txt file.  Exiting.')
                    return
            end
            
            motion_data{i,k}.frame_removal = cat(1,motion_data{i,k}.frame_removal, ...
                power_2014_motion_censoring(FD_vector,...
                                            DVAR_vector,...
                                            motion_data{i,k}.FD_threshold,...
                                            motion_data{i,k}.DVAR_threshold,...
                                            motion_data{i,k}.skip,...
                                            expected_contiguous_frame_count...
                )...
            );
            
            FD_all_runs = cat(1,FD_all_runs,FD_vector);
            DVAR_all_runs = cat(1,DVAR_all_runs,DVAR_vector);
        end
        
        motion_data{i,k}.format_string = format_generator(motion_data{i,k}.frame_removal);
        motion_data{i,k}.total_frame_count = length(motion_data{i,k}.frame_removal);
        motion_data{i,k}.remaining_frame_count = motion_data{i,k}.total_frame_count - sum(motion_data{i,k}.frame_removal);
        motion_data{i,k}.remaining_seconds = motion_data{i,k}.remaining_frame_count * motion_data{i,k}.epi_TR;
        motion_data{i,k}.remaining_frame_mean_FD = mean(FD_all_runs(motion_data{i,k}.frame_removal==0));
        motion_data{i,k}.remaining_frame_mean_DVAR = mean(DVAR_all_runs(motion_data{i,k}.frame_removal==0));
        
    end
end

save('power_2014_motion.mat','motion_data');
