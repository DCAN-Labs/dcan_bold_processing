%% Main Motion Censoring Function

function logical_removal_vector = power_2014_motion_censoring(FD_vector,DVAR_vector,FD_threshold,DVAR_threshold,skip,expected_contiguous_frame_count)

% Caluclate and combine the FD and DVAR removal vectors
FD_removal_vector = censor(FD_vector,FD_threshold);
DVAR_removal_vector = censor(DVAR_vector,DVAR_threshold);
prelim_removal_vector = logical(FD_removal_vector + DVAR_removal_vector);

logical_removal_vector = contiguous_frame_censor(prelim_removal_vector,expected_contiguous_frame_count);

% Skip cannot equal zero
if skip == 0
    skip = 1;
end

logical_removal_vector(1:skip) = true;

end

%% Censoring helper function

function logical_removal_vector = censor(motion_vector,motion_threshold)

logical_removal_vector = logical(motion_vector > motion_threshold);

end

%% Contiguous frame censoring helper function

function logical_removal_vector = contiguous_frame_censor(removal_vector,minimum_contiguous_frame_count)

contiguous_removal_vector = zeros(size(removal_vector));
contiguous_groups = bwlabel( logical((removal_vector - 1) * -1) );

for group = 1:max(contiguous_groups)
    if sum(contiguous_groups == group) < minimum_contiguous_frame_count
        contiguous_removal_vector(contiguous_groups == group) = true;
    else
        contiguous_removal_vector(contiguous_groups == group) = false;
    end
end

logical_removal_vector = logical(removal_vector + contiguous_removal_vector);

end