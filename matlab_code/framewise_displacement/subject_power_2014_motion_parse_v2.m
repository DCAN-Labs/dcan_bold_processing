function subject_power_2014_motion_parse_v2(movement_files,skip_seconds,epi_TR,expected_contiguous_frame_count)
% FD parsing Nipype wrapper pulling in list of motion files
% Modified to read DVAR_pre_reg, DVAR_post_reg, and DVAR_post_all and store
% values appropriately in the 51x51x3 cell array of structs
% 
% All variables (skips, TR, CurrentDirectory) set within Nipype
% Iterates over range of values from 0.1 to maximum of 0.5 FD threshold
% Outputs to .mat file
% FD_movement_files = {'motion_numbers0.txt', 'motion_numbers1.txt', 'motion_numbers2.txt'};
% FD_movement_files can be TXT or DAT files
% varargin (brain_radius_in_mm) only needs to be provided for DAT file inputs

FD_range = 0:0.01:0.5;
DVAR_range = 0:1:50;
motion_data = cell(length(FD_range),length(DVAR_range),3);
run_count = size(movement_files,2);

for i = 1:length(FD_range)
    for j = 1:length(DVAR_range)
        for k = 1:length({'DVAR Pre-Regression','DVAR Post-Regression','DVAR Post-All'})
            
            motion_data{i,j,k}.skip = ceil(skip_seconds/epi_TR); % Frames to skip
            motion_data{i,j,k}.epi_TR = epi_TR;   % EPI/BOLD/Functional Run Repetition Time
            motion_data{i,j,k}.FD_threshold = FD_range(i); % FD Threshold
            motion_data{i,j,k}.DVAR_threshold = DVAR_range(j); % DVAR Threshold
            motion_data{i,j,k}.frame_removal = []; % Initialize empty frame_removal vector
            
            FD_all_runs = [];
            DVAR_pre_reg_all_runs = [];
            DVAR_post_reg_all_runs = [];
            DVAR_post_all_all_runs = [];
            
            for run = 1:run_count
                movement_file  = movement_files{run};
                
                % Pad FD vector with a false frame (difference of frame "0" and 1)
                % before the first actual FD number (the difference of frame 1
                % and 2)
                switch movement_file(end-3:end)
                    case '.txt' % Read FD values from motion numbers files
                        fid = fopen(movement_file);
                        motion_cell = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f %f %f','HeaderLines',1);% Read FD values from motion numbers files
                        fclose(fid);
                        FD_vector = cat(1,0,motion_cell{1});
                        DVAR_pre_reg_vector = cat(1,0,motion_cell{2});
                        DVAR_post_reg_vector = cat(1,0,motion_cell{3});
                        DVAR_post_all_vector = cat(1,0,motion_cell{4});
                    otherwise
                        disp('Input is not a .txt file.  Exiting.')
                        return
                end
                
                switch k
                    case 1
                        motion_data{i,j,k}.frame_removal = cat(1,motion_data{i,j,k}.frame_removal, ...
                            power_2014_motion_censoring(FD_vector,...
                            DVAR_pre_reg_vector,...
                            motion_data{i,j,k}.FD_threshold,...
                            motion_data{i,j,k}.DVAR_threshold,...
                            motion_data{i,j,k}.skip,...
                            expected_contiguous_frame_count...
                            )...
                            );
                    case 2
                        motion_data{i,j,k}.frame_removal = cat(1,motion_data{i,j,k}.frame_removal, ...
                            power_2014_motion_censoring(FD_vector,...
                            DVAR_post_reg_vector,...
                            motion_data{i,j,k}.FD_threshold,...
                            motion_data{i,j,k}.DVAR_threshold,...
                            motion_data{i,j,k}.skip,...
                            expected_contiguous_frame_count...
                            )...
                            );
                    case 3
                        motion_data{i,j,k}.frame_removal = cat(1,motion_data{i,j,k}.frame_removal, ...
                            power_2014_motion_censoring(FD_vector,...
                            DVAR_post_all_vector,...
                            motion_data{i,j,k}.FD_threshold,...
                            motion_data{i,j,k}.DVAR_threshold,...
                            motion_data{i,j,k}.skip,...
                            expected_contiguous_frame_count...
                            )...
                            );
                end
                
                FD_all_runs = cat(1,FD_all_runs,FD_vector);
                DVAR_pre_reg_all_runs = cat(1,DVAR_pre_reg_all_runs,DVAR_pre_reg_vector);
                DVAR_post_reg_all_runs = cat(1,DVAR_post_reg_all_runs,DVAR_post_reg_vector);
                DVAR_post_all_all_runs = cat(1,DVAR_post_all_all_runs,DVAR_post_all_vector);
            end
            
            motion_data{i,j,k}.format_string = format_generator(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.total_frame_count = length(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.remaining_frame_count = motion_data{i,j,k}.total_frame_count - sum(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.remaining_seconds = motion_data{i,j,k}.remaining_frame_count * motion_data{i,j,k}.epi_TR;
            motion_data{i,j,k}.remaining_frame_mean_FD = mean(FD_all_runs(motion_data{i,j,k}.frame_removal==0));
            
            switch k
                case 1
                    motion_data{i,j,k}.remaining_frame_mean_DVAR = mean(DVAR_pre_reg_all_runs(motion_data{i,j,k}.frame_removal==0));
                case 2
                    motion_data{i,j,k}.remaining_frame_mean_DVAR = mean(DVAR_post_reg_all_runs(motion_data{i,j,k}.frame_removal==0));
                case 3
                    motion_data{i,j,k}.remaining_frame_mean_DVAR = mean(DVAR_post_all_all_runs(motion_data{i,j,k}.frame_removal==0));
            end
            
        end
    end
end

save('power_2014_motion.mat','motion_data');
