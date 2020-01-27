function subject_outliers_parse_BIDS(wb_command,motion_data_file,dtseries,output_dir,output_prefix)
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

load(motion_data_file);
FD_range_max = motion_data{length(motion_data)}.FD_threshold;
FD_range = 0:0.01:FD_range_max;

for i = 1:length(motion_data)

    motion_data{i}.fd_removal = motion_data{i}.frame_removal;
    motion_data{i}.outlier_removal = double(outlier_censoring(wb_command, dtseries, motion_data{i}.frame_removal));
    motion_data{i}.combined_removal = double(motion_data{i}.fd_removal | motion_data{i}.outlier_removal);

    motion_data{i}.fd_format_string = format_generator(motion_data{i}.fd_removal);
    motion_data{i}.outlier_format_string = format_generator(motion_data{i}.outlier_removal);
    motion_data{i}.combined_format_string = format_generator(motion_data{i}.combined_removal);

    motion_data{i}.remaining_fd_count = motion_data{i}.total_frame_count - sum(motion_data{i}.fd_removal);
    motion_data{i}.remaining_outlier_count = motion_data{i}.total_frame_count - sum(motion_data{i}.outlier_removal);
    motion_data{i}.remaining_combined_count = motion_data{i}.total_frame_count - sum(motion_data{i}.combined_removal);

    motion_data{i}.remaining_fd_seconds = motion_data{i}.remaining_fd_count * motion_data{i}.epi_TR;
    motion_data{i}.remaining_outlier_seconds = motion_data{i}.remaining_outlier_count * motion_data{i}.epi_TR;
    motion_data{i}.remaining_combined_seconds = motion_data{i}.remaining_combined_count * motion_data{i}.epi_TR;

end

save([output_dir filesep output_prefix '_outliers_power_2014_FD_only.mat'],'motion_data');
