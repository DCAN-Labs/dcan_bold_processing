function x = calc_FD_HCP_octave(prefix,brain_radius_in_mm)
if nargin < 2
    brain_radius_in_mm=50;
end

%% Read the first 6 columns

fileID = fopen(prefix,'r');

tline = fgetl(fileID);
X = [];
while ischar(tline)
	tline = strtrim(regexprep(tline,' +',' '));
    row = str2double(strsplit(tline, ' '));
    row = row(1:6);
    X=[X;row];
    tline = fgetl(fileID);
end

fclose(fileID);


%% Do the job
X(:,4:end)=X(:,4:end)*pi*brain_radius_in_mm/180; % Calculate length of arc in mm

dX=diff(X); % calculate derivatives
x=(sum(abs(dX),2)); % calculate FD


