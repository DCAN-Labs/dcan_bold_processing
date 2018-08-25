function value = read_power_2014_FD_only(motion_file,FD_threshold,key)

if FD_threshold > 0.5
    disp(['Input FD threshold (' num2str(FD_threshold) ') is greater than the maximum tested against (0.5)'])
    return
end

load(motion_file);
motion_struct = motion_data{floor(100*FD_threshold)+1};

if strcmp(key,'format_string')
    value = motion_struct.format_string;
elseif strcmp(key,'remaining_seconds')
    value = motion_struct.remaining_seconds;
elseif strcmp(key,'remaining_frame_count')
    value = motion_struct.remaining_frame_count;
elseif strcmp(key,'total_frame_count')
    value = motion_struct.total_frame_count;
elseif strcmp(key,'frame_removal')
    value = motion_struct.frame_removal;
elseif strcmp(key,'epi_TR')
    value = motion_struct.epi_TR;
elseif strcmp(key,'skip')
    value = motion_struct.skip;
else
    disp(['Key provided (' key ') does not match available keys.'])
end
