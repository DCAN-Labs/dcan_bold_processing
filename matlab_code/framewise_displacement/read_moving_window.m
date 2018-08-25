function value = read_moving_window(moving_window_file,window_index,key)


load(moving_window_file);
moving_window_struct = window_cell{window_index};

if strcmp(key,'format_string')
    value = moving_window_struct.format_string;
% elseif strcmp(key,'remaining_seconds')
%     value = moving_window_struct.remaining_seconds;
% elseif strcmp(key,'remaining_frame_count')
%     value = moving_window_struct.remaining_frame_count;
% elseif strcmp(key,'total_frame_count')
%     value = moving_window_struct.total_frame_count;
elseif strcmp(key,'frame_removal')
    value = moving_window_struct.frame_removal;
% elseif strcmp(key,'epi_TR')
%     value = moving_window_struct.epi_TR;
% elseif strcmp(key,'skip')
%     value = moving_window_struct.skip;
else
    disp(['Key provided (' key ') does not match available keys.'])
end
