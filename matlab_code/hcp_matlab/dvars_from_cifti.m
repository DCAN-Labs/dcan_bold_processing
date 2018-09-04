function DVAR=dvars_from_cifti(X);
% This  function calculates dvars based on grayordinates (WM and non brain excluded)

%% Check size
[g tr]=size(X);
if g<tr
    X=X';
end
dx=diff(X,[],2);
DVAR=sqrt(mean(dx.^2));
