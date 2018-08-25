function specific_filtered_movement_regressors(path_mov_reg,TR,filt_option,order,HCP_mat_path,FNL_preproc_version,PMU_dir,functype,restnum)
%% This function applies a low pas filter to the motion numbers and save it with a different number

%save('/mnt/max/home/andersp/FNLpreproc_v3_workspace.mat')
%load('/mnt/max/home/andersp/FNLpreproc_v3_workspace.mat')

% %% Input arguments
% path_mov_reg        Path until MNINonLinear/Results
% TR                  TR, in seconds
% filt_option         filter option, options 1:6, see below
% order               order of the filter to be applied
% LP_freq_min         Low pass frequency (in minutes) to be filtered
%% defaults
head_ratio_cm=5;
if exist('HCP_mat_path','var') == 0
    HCP_mat_path='/mnt/max/shared/utilities/HCP_Matlab';
end
%% paths
%addpath(genpath(HCP_mat_path));

%% when compiling comment out paths and uncomment this section
TR = str2num(TR);
filt_option = str2num(filt_option);
order = str2num(order);

% convert restnum to string to get movement_regressors path
restnum = num2str(restnum); % convert restnum to string for file path

%% Extract PMU data
%TODO: if PMU_dir does not equal 'false' then transformPMUOutput, extract RR frequency, and re-calculate fc_RR_min and fc_RR_max used to create the notch filter.
TransformPMUOutput(PMU_dir);
[RRa_Hz, HRa_Hz, RR_Hz, HR_Hz]=read_aliased_PMU_data(PMU_dir, TR);
% calculate and overwrite fc_RR_min from RR_min here
RR_min=RR_Hz*60;
HR_min=HR_Hz*60;
fc_RR_min=RR_min-4;
fc_RR_max=RR_min+4;
disp(['fc_RR_min = ' num2str(fc_RR_min)]);
disp(['fc_RR_max = ' num2str(fc_RR_max)]);

%% filter design

filt_method{1}='None';
filt_method{2}='None_and_ignore_Y';
filt_method{3}='Filt_all';
filt_method{4}='Filt_only_Y';
filt_method{5}='FiltFilt_all';
filt_method{6}='FiltFilt_only_Y';

%handle niftis that store TR in ms
if (TR > 20)
    TR=TR/1000;
end

% Create the notch filter
fc_RR_bw=[fc_RR_min,fc_RR_max];
rr=fc_RR_bw/60;

fs = 1/TR;
fNy=fs/2;

fa=abs(rr-floor((rr+fNy)/fs)*fs);

W_notch = fa/fNy;    
Wn=mean(W_notch);
bw=abs(diff(W_notch));
[b_filt,a_filt]=iirnotch(Wn,bw);
num_f_apply=floor(order/2); % if order<4 apply filter 1x, if order=4 2x, if order=6 3x


%% Read individual movement regressors files
f=filesep;
%defined_path=strtrim(ls([path_mov_reg f 'REST*' f 'Movement_Regressors.txt']));
%defined_path=regexp(defined_path,'\n','split'); %-- FOR MATLAB

% sorted_path=sort_RS(defined_path);
%n=size(defined_path,2);



% Read motion numbers
file_mov_reg=[path_mov_reg f functype restnum f 'Movement_Regressors.txt' ]; %-- for MATLAB USE
disp(file_mov_reg)
%file_mov_reg=defined_path(i,1:find(defined_path(i,:)=='.')+3);
MR = dlmread(file_mov_reg);
MR_ld=make_friston_regressors(MR);%% Using this function to only get the linear displacements
MR_ld=MR_ld(:,1:6);


switch filt_option
    case 1 %'None';
        MR_filt=MR_ld;
            
    case 2 %'None_and_ignore_Y';
        MR_filt=MR_ld;
        MR_filt(:,2)=0;
            
    case 3 %'Filt_all';
        MR_filt=filter(b_filt,a_filt,MR_ld);
        for i=1:num_f_apply-1
            MR_filt=filter(b_filt,a_filt,MR_filt);
        end
            
    case 4 %'Filt_only_Y';
        MR_filt=MR_ld;
        MR_filt(:,2)=filter(b_filt,a_filt,MR_ld(:,2));
        for i=1:num_f_apply-1
            MR_filt(:,2)=filter(b_filt,a_filt,MR_filt(:,2));
        end
            
    case 5 %'FiltFilt_all';
        MR_filt=filtfilt(b_filt,a_filt,MR_ld);
        for i=1:num_f_apply-1
            MR_filt=filtfilt(b_filt,a_filt,MR_filt);
        end
            
    case 6 %'FiltFilt_only_Y';
        MR_filt=MR_ld;
        MR_filt(:,2)=filtfilt(b_filt,a_filt,MR_ld(:,2));
        for i=1:num_f_apply-1
            MR_filt(:,2)=filtfilt(b_filt,a_filt,MR_filt(:,2));
        end    
end


hd_mm=10*head_ratio_cm; % cm to mm conversion
MR_backed=MR_filt;
MR_backed(:,4:end)=180*MR_backed(:,4:end)/(pi*hd_mm);
second_derivs=cat(1,[0 0 0 0 0 0], diff(MR_backed(:,1:6)));


%%
order_string = num2str(order);
if strcmp(order_string,'')
order_string = '[]';
end
[path_orig file_orig ext_orig]=fileparts(file_mov_reg);
FNL_preproc_version=char(FNL_preproc_version)
file_filtered=[file_orig '_' FNL_preproc_version];
full_file_filtered=[path_orig f file_filtered ext_orig];
disp(full_file_filtered)
dlmwrite(full_file_filtered,MR_backed,' ');



    %NOTES: THIS has been changed to dlmwrite to work with octave. Noted by
%Eric Feczko on 3/8/17 at 5:01 p.m.
%    fileID = fopen(full_file_filtered,'w');
%    
%   
%    num_fmt='%16.6f';
%    formatSpec = [num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt '\n'];
%    for j=1:size(MR_backed,1)
%        fprintf(fileID,formatSpec, MR_backed(j,1), MR_backed(j,2), MR_backed(j,3), MR_backed(j,4), MR_backed(j,5), MR_backed(j,6));
%    end
%    fclose(fileID);
    

%function sorted_path=sort_RS(defined_path)
%aa=dir(defined_path);
%
%nn=size(aa,1);
%ix_ls=zeros(nn,1);
%sorted_path=cell(nn,1);
%for ii=1:nn
%    foo=strsplit(aa(ii).folder,filesep);
%    try
%        bb=foo{end};
%        ix_ls(ii)=str2num(bb(5:end));
%    catch
%        bb=foo{end-1};
%        ix_ls(ii)=str2num(bb(5:end));
%    end
%end
%[cc, ix]=sort(ix_ls);
%
%for ii=1:nn
%    aa(ix(ii)).folder;
%    sorted_path{ii}=[aa(ix(ii)).folder filesep aa(ix(ii)).name];
%end

exit
