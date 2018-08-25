function motion = motionfile_parse(prefix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION: SPI April 10th 2014

% prefix is the directory and file name of the movement file for each run
% EG: W:\NIGG\Analyses\Pilot\processed_data\s0053-1/movement/s0053-1_b1_faln_dbnd_xr3d.dat
% x is the FD for a given run without removing the skips
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Finding the index for the footer text
clear motion1
fid=fopen(prefix);
motion1=textscan(fid,'%s%s%s%s%s%s%s%s','HeaderLines',3);
val=strcmp(motion1{1,1},'#counting');
fclose(fid);

M=find(val==1);
N=M-1;

%loading the dat file into matlab

fid=fopen(prefix);
motion=textscan(fid,'%d%f%f%f%f%f%f%f', N,'HeaderLines',3);
fclose(fid);

