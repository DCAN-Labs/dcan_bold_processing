function value = read_FD(FD_file,FD_threshold,key)

if FD_threshold > 0.5
    disp(['Input FD threshold (' num2str(FD_threshold) ') is greater than the maximum tested against (0.5)'])
    return
end

load(FD_file);
FD_struct = FD_data{floor(100*FD_threshold)+1};

if strcmp(key,'format_string')
    value = FD_struct.format_string;
elseif strcmp(key,'remaining_seconds')
    value = FD_struct.remaining_seconds;
elseif strcmp(key,'remaining_frame_count')
    value = FD_struct.remaining_frame_count;
elseif strcmp(key,'total_frame_count')
    value = FD_struct.total_frame_count;
elseif strcmp(key,'frame_removal')
    value = FD_struct.frame_removal;
elseif strcmp(key,'epi_TR')
    value = FD_struct.epi_TR;
elseif strcmp(key,'skip')
    value = FD_struct.skip;
else
    disp(['Key provided (' key ') does not match available keys.'])
end
