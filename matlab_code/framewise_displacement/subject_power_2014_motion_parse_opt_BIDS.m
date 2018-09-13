function subject_power_2014_motion_parse_opt_BIDS(movement_files,skip_seconds,epi_TR,expected_contiguous_frame_count,output_dir,output_prefix)
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

N_FD=51;
N_DVAR=51;

max_FD=0.5;
max_DVAR=50;

% FD_range = str2num(num2str(linspace(0,max_FD,N_FD)));
% DVAR_range = str2num(num2str(linspace(0,max_DVAR,N_DVAR)));

FD_range = linspace(0,max_FD,N_FD);
DVAR_range = linspace(0,max_DVAR,N_DVAR);

clear motion_data

motion_data = cell(N_FD,N_DVAR,3);

run_count = length(movement_files);
FD_all_runs = [];
DVAR_pre_reg_all_runs = [];
DVAR_post_reg_all_runs = [];
DVAR_post_all_all_runs = [];
ix=[];
for run = 1:run_count
    movement_file  = movement_files{run};
    
    % Pad FD vector with a false frame (difference of frame "0" and 1)
    % before the first actual FD number (the difference of frame 1
    % and 2)
    clear motion_cell
    switch movement_file(end-3:end)
        case '.txt' % Read FD values from motion numbers files
            fid = fopen(movement_file);
            motion_cell = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f %f %f','HeaderLines',1);% Read FD values from motion numbers files
            fclose(fid);
            FD_vector = cat(1,0,motion_cell{1});
            DVAR_pre_reg_vector = cat(1,0,motion_cell{2});
            DVAR_post_reg_vector = cat(1,0,motion_cell{3});
            DVAR_post_all_vector = cat(1,0,motion_cell{4});
            
            FD_all_runs = cat(1,FD_all_runs,FD_vector);
            DVAR_pre_reg_all_runs = cat(1,DVAR_pre_reg_all_runs,DVAR_pre_reg_vector);
            DVAR_post_reg_all_runs = cat(1,DVAR_post_reg_all_runs,DVAR_post_reg_vector);
            DVAR_post_all_all_runs = cat(1,DVAR_post_all_all_runs,DVAR_post_all_vector);
            ix=cat(1,ix,repmat(run,size(FD_vector)));
        otherwise
            disp('Input is not a .txt file.  Exiting.')
            return
    end
end
DVARS=[DVAR_pre_reg_all_runs DVAR_post_reg_all_runs DVAR_post_all_all_runs];

for i = 1:N_FD
    for j = 1:N_DVAR
        for k = 1:length({'DVAR Pre-Regression','DVAR Post-Regression','DVAR Post-All'})
            
            motion_data{i,j,k}.skip = ceil(skip_seconds/epi_TR); % Frames to skip
            motion_data{i,j,k}.epi_TR = epi_TR;   % EPI/BOLD/Functional Run Repetition Time
            motion_data{i,j,k}.FD_threshold = FD_range(i); % FD Threshold
            motion_data{i,j,k}.DVAR_threshold = DVAR_range(j); % DVAR Threshold
            motion_data{i,j,k}.frame_removal = []; % Initialize empty frame_removal vector
            
            for ii=1:run
                motion_data{i,j,k}.frame_removal = cat(1,motion_data{i,j,k}.frame_removal, ...
                power_2014_motion_censoring(FD_all_runs(ix==ii),...
                DVARS(ix==ii,k),...
                motion_data{i,j,1}.FD_threshold,...
                motion_data{i,j,1}.DVAR_threshold,...
                motion_data{i,j,1}.skip,...
                expected_contiguous_frame_count...
                )...
                );
            end
            
            motion_data{i,j,k}.format_string = format_generator(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.total_frame_count = length(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.remaining_frame_count = motion_data{i,j,k}.total_frame_count - sum(motion_data{i,j,k}.frame_removal);
            motion_data{i,j,k}.remaining_seconds = motion_data{i,j,k}.remaining_frame_count * motion_data{i,j,k}.epi_TR;
            motion_data{i,j,k}.remaining_frame_mean_FD = mean(FD_all_runs(motion_data{i,j,k}.frame_removal==0));
            motion_data{i,j,k}.remaining_frame_mean_DVAR = mean(DVARS(motion_data{i,j,k}.frame_removal==0,k));
        end
    end
end

save([output_dir filesep output_prefix '_power_2014_motion.mat'],'motion_data');
