function R=conc_and_detrend(wm, vent, WB, FR)
%% function R=conc_and_detrend(wm, vent, WB, FR)
%
% This function concatenate and detrend regressors. It adds derivatives to
% wm, vent and WB
%% 
sR=[wm vent WB]; %simple regressors
dR=diff(sR);
dR=[0 0 0;dR];
R=[sR dR FR];
R=detrend(R);