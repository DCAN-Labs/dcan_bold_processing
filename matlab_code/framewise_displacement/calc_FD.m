function x = calc_FD(prefix,brain_radius_in_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION: SPI April 10th 2014, commented by Eric Earl 9/8/2015
% 
% INPUTS:
% prefix is the directory and file name of the movement file per run
%  - e.g. W:/NIGG/Pilot/processed_data/subject1/movement/subject1_b1_faln_dbnd_xr3d.dat
% brain_radius_in_mm should be the population average brain radius in
% millimeters
%  - e.g. 50 (for human youth or adult)
% OUTPUTS:
% x is the FD for a given run without removing the skips
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Finding the index for the footer text
clear motion1 % starting fresh
fid=fopen(prefix); % open the file
motion1=textscan(fid,'%s%s%s%s%s%s%s%s','HeaderLines',3); % read in the 8 columns skipping the 3 header lines
val=strcmp(motion1{1,1},'#counting'); % look for the line with "#counting" on it
fclose(fid); % close the file

M=find(val==1); % find the line where there is "#counting" in the text
N=M-1; % N=number of frames is M - 1

%loading the dat file into matlab

fid=fopen(prefix); % open the file again
motion=textscan(fid,'%d%f%f%f%f%f%f%f', N,'HeaderLines',3); % this time parse the file correctly with floats
fclose(fid); % close the file

%% OMD: Consider using this code:
% X=zeros(N,6);
% for i=1:6
%     X(:,i)=motion{i+1};
% end
% X(:,4:end)=X(:,4:end)*pi*brain_radius_in_mm/180;
% 
% dX=diff(X);
% x=(sum(abs(dX),2));
%%

%Converting the rotation into mm
% rotation is the result of converting from degrees to radians by
% multiplying by 2*pi radians and dividing by 360 degrees
x_rot=(brain_radius_in_mm*pi/180)*(motion{1,5});
y_rot=(brain_radius_in_mm*pi/180)*(motion{1,6});
z_rot=(brain_radius_in_mm*pi/180)*(motion{1,7});

% clear out any leftovers from last time
clear diff_x
clear diff_y
clear diff_z
clear rot_diff_x
clear rot_diff_y
clear rot_diff_z

dif=zeros(6,N-1); % preallocate the difference matrix
for b=2:N % do the frame-to-frame differences
    dif(1,b-1)=motion{1,2}(b-1)-motion{1,2}(b);
    dif(2,b-1)=motion{1,3}(b-1)-motion{1,3}(b);
    dif(3,b-1)=motion{1,4}(b-1)-motion{1,4}(b);
    dif(4,b-1)=x_rot(b-1)-x_rot(b);
    dif(5,b-1)=y_rot(b-1)-y_rot(b);
    dif(6,b-1)=z_rot(b-1)-z_rot(b);
end

x=(sum(abs(dif),1))'; % return the transpose of the sum of the absolute values of all the differences frame pair by frame pair
