function motion_summary(motion_numbers_file,FD_thresholds_file,output_folder)

load(motion_numbers_file)
load(FD_thresholds_file)

motion_indices = 6:5:51;
remaining_seconds = zeros(1,length(motion_indices));
FD_thresholds = zeros(1,length(motion_indices));

FD = motion_numbers.FD;
DVAR_pre_reg = motion_numbers.DVAR_pre_reg;
DVAR_post_reg = motion_numbers.DVAR_post_reg;
DVAR_post_all = motion_numbers.DVAR_post_all;
keep_frames = FD~=0;

outstruct.FD_mean = mean(FD(keep_frames));
outstruct.DVAR_pre_regression_mean = mean(DVAR_pre_reg(keep_frames));
outstruct.DVAR_post_regression_mean = mean(DVAR_post_reg(keep_frames));
outstruct.DVAR_post_all_mean = mean(DVAR_post_all(keep_frames));

outstruct.total_time_in_seconds = motion_data{1}.epi_TR*motion_data{1}.total_frame_count;

for i = 1:length(motion_indices)
    remaining_seconds(i) = motion_data{motion_indices(i)}.remaining_seconds;
    FD_thresholds(i) = motion_data{motion_indices(i)}.FD_threshold;
end

for i = 1:length(motion_indices)
    outstruct.(strrep([ 'remaining_seconds_at_FD_' num2str(FD_thresholds(i),'%.2f') ],'.','_')) = remaining_seconds(i);
end

writetable(struct2table(outstruct),[output_folder filesep 'motion_summary_table.csv']);
