function R=conc(wm, vent, FR)
%% function R=conc(wm, vent, FR)
%
% This function concatenate and detrend regressors. It adds derivatives to
% wm, vent
%% 
sR=[wm vent]; %simple regressors
dR=diff(sR);
dR=[0 0;dR];
R=[sR dR FR];
% R=detrend(R);
