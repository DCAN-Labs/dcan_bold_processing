function filtered_movement_regressors(path_mov_reg, TR, filt_option, order, LP_freq_min, filt_type, fc_RR_min, fc_RR_max, output)
%% This function applies a low pas filter to the motion numbers and save it with a different number

% %% Input arguments
% path_mov_reg        Path until MNINonLinear/Results
% TR                  TR, in seconds
% filt_option         filter option, options 1:6, see below
% order               order of the filter to be applied
% LP_freq_min         Low pass frequency (in minutes) to be filtered
%% defaults
head_ratio_cm = 5;

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

switch filt_type

    case 'lp'
        
        hr_min = LP_freq_min; % recasted to reuse code
        hr = hr_min/60;

        fs = 1/TR;
        fNy = fs/2;

        fa = abs(hr - floor((hr + fNy) / fs) * fs);

        % cutting frequency normalized between 0 and nyquist
        Wn = min(fa)/fNy;
        if ~isempty(order)
            b_filt = fir1(order, Wn, 'low');
            a_filt = 1;
        end
        num_f_apply = 0;

    case 'notch'
        
        fc_RR_bw = [fc_RR_min, fc_RR_max];
        rr = fc_RR_bw / 60;

        fs = 1 / TR;
        fNy = fs / 2;

        fa = abs(rr - floor((rr+fNy) / fs) * fs);

        W_notch = fa / fNy;
        Wn = mean(W_notch);
        bw = diff(W_notch);
        [b_filt, a_filt] = iirnotch(Wn, bw);
        num_f_apply = floor(order / 2); % if order<4 apply filter 1x, if order=4 2x, if order=6 3x

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
            for i=1:num_f_apply-1
                MR_filt = filter(b_filt,a_filt,MR_filt);
            end
            
        case 4 %'Filt_only_Y';
            MR_filt = MR_ld;
            MR_filt(:,2) = filter(b_filt,a_filt,MR_ld(:,2));
            for i=1:num_f_apply-1
                MR_filt(:,2) = filter(b_filt,a_filt,MR_filt(:,2));
            end
            
        case 5 %'FiltFilt_all';
            MR_filt = filtfilt(b_filt,a_filt,MR_ld);
            for i=1:num_f_apply-1
                MR_filt = filtfilt(b_filt,a_filt,MR_filt);
            end
            
        case 6 %'FiltFilt_only_Y';
            MR_filt = MR_ld;
            MR_filt(:,2) = filtfilt(b_filt,a_filt,MR_ld(:,2));
            for i=1:num_f_apply-1
                MR_filt(:,2) = filtfilt(b_filt,a_filt,MR_filt(:,2));
            end
    end
    
    
    hd_mm = 10 * head_ratio_cm; % cm to mm conversion
    MR_backed = MR_filt;
    MR_backed(:, 4:end) = 180*MR_backed(:,4:end)/(pi*hd_mm);

    %% The last 6 movement regressors are derivatives of the first 6. Needed for task fMRI, but not used in the Fair Lab
    second_derivs=cat(1, [0 0 0 0 0 0], diff(MR_backed(:,1:6)));
    MR_backed = [MR_backed second_derivs];
   
    % [path_orig, file_orig, ext_orig] = fileparts(file_mov_reg);
    dlmwrite(output, MR_backed, ' ');
    
end

exit
