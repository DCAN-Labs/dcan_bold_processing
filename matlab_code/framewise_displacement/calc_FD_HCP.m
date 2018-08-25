function x = calc_FD_HCP(prefix,brain_radius_in_mm)
if nargin < 2
    brain_radius_in_mm=50;
end

%% Read the first 6 columns

formatSpec = '%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(prefix,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
X = [dataArray{1:end-1}];
%% Do the job

X(:,4:end)=X(:,4:end)*pi*brain_radius_in_mm/180; % Calculate length of arc in mm

dX=diff(X); % calculate derivatives
x=(sum(abs(dX),2)); % calculate FD

