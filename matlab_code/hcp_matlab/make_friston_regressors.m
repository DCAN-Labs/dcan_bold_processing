function FR=make_friston_regressors(MR, head_ratio_cm)

%% function R=make_friston_regressors(R, head_ratio_mm)
% This function takes as input a matrix having 6 columns, corresponding to
% the 6DOF movement correction parameters of EPI data and then it
% calculates the 24 corresponding Friston regressors, as reported in
%
%% Assumptions
%
% First 3 columns of R are in mm. Columns 4-6 are in degrees, as calculated
% by the HCP pipelines
%% Mandatory Inputs
%
% MR
%    A matrix of size (r x c), where r is the number of time points and
%    c are the 6 DOF regresors. If the number of columns is more than 6, the
%    code only considers the first 6 columns
%% Optional inputs
%
% head_ratio_cm
%   Head ratio, in cm. If not provided, the code use 5 cm as default
%% Outputs
%
%   FR
%       24 movement regresors, which corresponds to the following operations:
%
%   FR=[MR MR.^2 MR_(t-1) MR_(t-1).^2] % Note, no detrending here
%% Temporary variables and check inputs

if nargin<2
    head_ratio_cm=5;
end
hd_mm=10*head_ratio_cm; % cm to mm conversion

MR=MR(:,1:6); % removing derivatives
%% calculating arc length

MR(:,4:end)=MR(:,4:end)*pi*hd_mm/180;
%% 
FR=[MR MR.^2];
dummy=FR;
dummy(2:end,:)=FR(1:end-1,:);
dummy(1,:)=0;
FR=[FR dummy];
