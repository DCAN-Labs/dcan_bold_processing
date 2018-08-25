function subject_FD_parse(FD_movement_files,skip_seconds,epi_TR,brain_radius_in_mm)
% FD parsing Nipype wrapper pulling in list of motion files
% All variables (skips, TR, CurrentDirectory) set within Nipype
% Iterates over range of values from 0.1 to maximum of 0.5 FD threshold
% Outputs to .mat file
% FD_movement_files = {'motion_numbers0.txt', 'motion_numbers1.txt', 'motion_numbers2.txt'};
% FD_movement_files can be TXT or DAT files
% varargin (brain_radius_in_mm) only needs to be provided for DAT file inputs
% NOTE: need to comment on how FD_movement_file is structured

FD_range = 0:0.01:0.5; 
FD_data = {};
run_count = size(FD_movement_files,2);

for i = 1:length(FD_range)
    
    FD_data{i}.skip = ceil(skip_seconds/epi_TR); % Frames to skip
    FD_data{i}.epi_TR = epi_TR;   % EPI/BOLD/Functional Run Repetition Time
    FD_data{i}.FD_threshold = FD_range(i); % FD Threshold
    FD_data{i}.frame_removal = []; % Initialize empty frame_removal vector
    
    FD_all_runs = [];
    
    for j = 1:run_count
        movement_file  = FD_movement_files{j};
        
        % Pad FD vector with a false frame (difference of frame "0" and 1) 
        % before the first actual FD number (the difference of frame 1
        % and 2)
        switch movement_file(end-3:end)
            case '.txt' % Read FD values from motion numbers files
                fid = fopen(movement_file);
                FD_cell = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*[^\n]','HeaderLines',1);% Read FD values from motion numbers files
                fclose(fid);
                FD_vector = cat(1,0,FD_cell{1});
            case '.dat' % Caluclate FD values from DAT files
                FD_vector = cat(1,0,calc_FD(movement_file,brain_radius_in_mm));
            otherwise
                disp('Input is neither a .dat or .txt file.  Exiting.')
                return
        end
        
        FD_data{i}.frame_removal = cat(1,FD_data{i}.frame_removal, ...
            pre_post_motion_censoring(FD_vector,FD_data{i}.FD_threshold,FD_data{i}.skip));
        
        FD_all_runs = cat(1,FD_all_runs,FD_vector);
        
    end
    
    FD_data{i}.format_string = format_generator(FD_data{i}.frame_removal);
    FD_data{i}.total_frame_count = length(FD_data{i}.frame_removal);
    FD_data{i}.remaining_frame_count = FD_data{i}.total_frame_count - sum(FD_data{i}.frame_removal);
    FD_data{i}.remaining_seconds = FD_data{i}.remaining_frame_count * FD_data{i}.epi_TR;
    FD_data{i}.remaining_frame_mean_FD = mean(FD_all_runs(FD_data{i}.frame_removal==0));
    
end

save('FD.mat','FD_data');
