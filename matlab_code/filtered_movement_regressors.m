function filtered_movement_regressors(path_mov_reg, TR, filt_option, order, LP_freq_min, filt_type, fc_RR_min, fc_RR_max, output)
%% This function applies a low pas filter to the motion numbers and save it with a different number

% %% Input arguments
% path_mov_reg        Path until MNINonLinear/Results
% TR                  TR, in seconds
% filt_option         filter option, options 1:6, see below
% order               order of the filter to be applied
% LP_freq_min         Low pass frequency (in minutes) to be filtered

%% when compiling comment out paths and uncomment this section
TR = str2num(TR);
filt_option = str2num(filt_option);
order = str2num(order);
LP_freq_min = str2num(LP_freq_min);
fc_RR_min = str2num(fc_RR_min);
fc_RR_max = str2num(fc_RR_max);

%% filter design

filt_method{1} = 'None';
filt_method{2} = 'None_and_ignore_Y';
filt_method{3} = 'Filt_all';
filt_method{4} = 'Filt_only_Y';
filt_method{5} = 'FiltFilt_all';
filt_method{6} = 'FiltFilt_only_Y';

if filt_option>4
    order = order/2; % because filtfilt does forward and backward filtering
end

fs = 1 / TR;
fNy = fs / 2;
switch filt_type
    case 'lp'
        hr_min = LP_freq_min; % recasted to reuse code
        hr = hr_min/60;
        [b_filt,a_filt] = butter(order,rr/fNy,'stop');[b_filt,a_filt] = butter(order/2,hr/fNy,'low');
    case 'notch'
        fc_RR_bw = [fc_RR_min, fc_RR_max];
        rr = fc_RR_bw / 60;
        [b_filt,a_filt] = butter(order,rr/fNy,'stop');
end


%% Read individual movement regressors files
if exist(path_mov_reg) == 7
    pathstring = [path_mov_reg filesep '*' filesep 'Movement_Regressors.txt'];
elseif exist(path_mov_reg) == 2
    pathstring = path_mov_reg;
else
    disp([path_mov_reg ' does not exist. Exiting...'])
    exit
end
path_contents = dir(pathstring);
n = size(path_contents,1);
for i=1:n
    
    % Read motion numbers
    file_mov_reg = [path_contents(i).folder filesep path_contents(i).name]; %-- for MATLAB USE
    MR = dlmread(file_mov_reg);
    MR_ld=make_friston_regressors(MR);%% Using this function to only get the linear displacements
    MR_ld=MR_ld(:,1:6);
    
    switch filt_option      
        case 1 %'None';
            MR_filt = MR_ld;
            
        case 2 %'None_and_ignore_Y';
            MR_filt = MR_ld;
            MR_filt(:,2) = 0;
            
        case 3 %'Filt_all';
            MR_filt = filter(b_filt,a_filt,MR_ld);
            
        case 4 %'Filt_only_Y';
            MR_filt = MR_ld;
            MR_filt(:,2) = filter(b_filt,a_filt,MR_ld(:,2));
            
        case 5 %'FiltFilt_all';
            MR_filt = filtfilt(b_filt,a_filt,MR_ld);
            for i=1:num_f_apply-1
                MR_filt = filtfilt(b_filt,a_filt,MR_filt);
            end
            
        case 6 %'FiltFilt_only_Y';
            MR_filt = MR_ld;
            MR_filt(:,2) = filtfilt(b_filt,a_filt,MR_ld(:,2));
    end
    
    MR_backed = MR_filt;

    %% The last 6 movement regressors are derivatives of the first 6. Needed for task fMRI, but not used in the Fair Lab
    second_derivs=cat(1, [0 0 0 0 0 0], diff(MR_backed(:,1:6)));
    MR_backed = [MR_backed second_derivs];
   
    % [path_orig, file_orig, ext_orig] = fileparts(file_mov_reg);
    dlmwrite(output, MR_backed, ' ');
    
end

exit
