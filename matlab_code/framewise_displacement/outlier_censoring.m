function output_removal_vector = outlier_censoring(wb_command, dtseries, input_removal_vector)

outlier_removal_vector = zeros(size(input_removal_vector,1), size(input_removal_vector,2));

[status,cmdout] = system([wb_command ' -cifti-stats ' dtseries ' -reduce STDEV']);

if status ~= 0
    return
end

% parse standard deviation values
stdevs = cell2mat(textscan(cmdout, '%f\n'));

% find the kept frames
input_idx = find(~input_removal_vector);

% find outliers
outlier_vector = isoutlier(stdevs(input_idx),'median');

% find indices of outliers within the kept frames
outlier_idx = find(outlier_vector);

% put the outlier flags into the outliers vector as ones (to censor)
outlier_removal_vector( input_idx(outlier_idx) ) = 1;

% convert outliers flags vector to a strict logical vector
output_removal_vector = logical(outlier_removal_vector);
