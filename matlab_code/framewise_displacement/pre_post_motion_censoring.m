function logical_removal_vector = pre_post_motion_censoring(FD_vector,FD_threshold,skip)
% Code from OHSU Fair Neuroimaging Lab
% Written (2014) and commented (9/8/2015) by Eric Earl
% 
% takes in a vector of FD values, an FD threshold, and a number of frames
% to skip and outputs a logical frame removal vector of 1s and 0s

% logical vector of frames to remove
FD_remove = FD_vector > FD_threshold;

% Remove adjacent frames using circshift
FD_remove_bkwd1 = circshift(FD_remove,-1);
FD_remove_fwd1 = circshift(FD_remove,1);
FD_remove_fwd2 = circshift(FD_remove,2);

% Fix the circular shifts added at the beginnings and ends
FD_remove_bkwd1(end) = 0;
FD_remove_fwd1(1) = 0;
FD_remove_fwd2(1:2) = 0;

% make a vector containing numbers greater than zero for frames to remove
removal_vector = sum(cat(2,FD_remove,FD_remove_bkwd1,FD_remove_fwd1,FD_remove_fwd2),2);

% Skip cannot equal zero for FD
if skip == 0
    skip = 1;
end

% insert the skip at the beginning
removal_vector(1:skip) = 1;

% finally output the logical vector of removals
logical_removal_vector = removal_vector > 0;
