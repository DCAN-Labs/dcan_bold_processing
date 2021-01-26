function dcan_signal_processing(config_path)
% DCAN Signal Processsing
%
%  runs motion filtering and signal regression on HCP grayordinate results.
%  accepts a configuration file (.json) which specifies parameters as well
%  as expected file paths.
%

%% load variables from config
conf_str = fileread(config_path);
conf_str = strrep(conf_str, ' ', ''); % hack to get loadjson to work
conf_json = loadjson(conf_str);

path_wb_c = conf_json.path_wb_c; %path to wb_command
bp_order= conf_json.bp_order;
lp_Hz= conf_json.lp_Hz;
hp_Hz= conf_json.hp_Hz;
TR= conf_json.TR;
fd_th= conf_json.fd_th;
path_cii= conf_json.path_cii;
path_ex_sum= conf_json.path_ex_sum;
FNL_preproc_CIFTI_basename= conf_json.FNL_preproc_CIFTI_basename;
fMRIName= conf_json.fMRIName;
file_wm= conf_json.file_wm;
file_vent= conf_json.file_vent;
file_mov_reg= conf_json.file_mov_reg;
motion_filename= conf_json.motion_filename;
skip_seconds= conf_json.skip_seconds;
result_dir= conf_json.result_dir;
int_method= 'linear';
%silence warnings
warning('off', 'all')

%% Predefining variables

%handle niftis that store TR in ms
if (TR > 20)
    TR=TR/1000;
end
disp(['TR = ' num2str(TR)]);
skip_frames = floor(skip_seconds/TR);

%% Read cifti
cifti_txt_path = [result_dir filesep 'temp_FNL_cifti.txt'];
cmd = [path_wb_c ' -cifti-convert -to-text ' path_cii ' ' cifti_txt_path];
disp(cmd)
system(cmd);
X = dlmread(cifti_txt_path);
system(['rm -f ' cifti_txt_path]);
DVAR_pre_reg=dvars_from_cifti(X);

%% readwm
fileid = fopen(file_wm);
wm = textscan(fileid, '%f');
wm = wm{1};
fclose(fileid);

%% read vent
fileid = fopen(file_vent);
vent = textscan(fileid, '%f');
vent = vent{1};
fclose(fileid);

%% read mov_reg and calculate Friston regressors
MR = dlmread(file_mov_reg);
MR = MR(:, [1 2 3 4 5 6]);
FR = make_friston_regressors(MR);

%% Calculate residuals (regression)
[r, c] = size(X');
Rr = zeros(r,c);

%% Create FD.mat
FD = calc_FD_HCP(file_mov_reg);
save([result_dir filesep 'FD.mat'], 'FD', '-v7')
txt_FD_file=[path_ex_sum filesep 'FD_' fMRIName '.txt'];
save(txt_FD_file, 'FD','-ascii')
th = [0; FD]<=fd_th;
th(1:skip_frames) = 0;
th = th==1;

%% Concatenate reg and detrened
WB = mean(X)';
try
    % 	    R=conc_and_detrend(wm, vent, WB, FR);
    
    R = conc(wm, vent, WB, FR);
    mR = mean(R(th, :));
    R = R-mR;
    n_reg = size(R, 2);
    bR = zeros(n_reg, 2);
    x = 1:r;
    local_x = x(th);
    for i=1:n_reg
        bR(i,:) = polyfit(local_x,R(th,i)', 1);
        R(:,i) = R(:,i) - polyval(bR(i,:), x)';
    end
    
    
catch ME
    disp(ME.message)
    disp(['Check your Movement_Regressors.txt for inconsistencies.  ' ...
          'Exiting...'])
    exit
end


%% detrend bold signal
mX = mean(X(:,th),2);
Xd = X - mX;
bXd = zeros(c,2);
for i=1:c
    bXd(i,:) = polyfit(local_x,Xd(i,th),1);
    Xd(i,:) = Xd(i,:) - polyval(bXd(i,:),x);
end
Xd = Xd';
w = 'stats:regress:RankDefDesignMat';
warning('off', w)


ts1 = Xd(th,:);
ts2 = R(th,:);

for i=1:c
    b = regress(ts1(:,i),ts2);
    Rr(:,i) = Xd(:,i) - R*b;
end
DVAR_post_reg = dvars_from_cifti(Rr);

%% remake cifti for residuals from regressors
cifti_resid_txt_path = [result_dir filesep 'temp_FNL_cifti_resid.txt'];

FNL_preproc_CIFTI_resid_name = strcat(FNL_preproc_CIFTI_basename, ...
                                      '_resid.dtseries.nii');
dlmwrite(cifti_resid_txt_path, Rr','delimiter' , ' ');

cmd = [path_wb_c ' -cifti-convert -from-text ' cifti_resid_txt_path ' ' ...
       path_cii ' ' result_dir filesep FNL_preproc_CIFTI_resid_name];

disp(cmd)
system(cmd);
system(['rm -f ' cifti_resid_txt_path]);

%% interpolation
ix_ol = find(th==0); % find index movers > fd_th
ix_in = find(th==1); % find index low motion number

if numel(ix_in) <= 10
    int_method = 'none';
    err_msg = ['WARNING: There are only ', num2str(numel(ix_in)), ...
        ' frames below the fd threshold of ', ...
        num2str(fd_th), '. Skipping interpolation.'];
    disp(err_msg)
end

Rr_int = Rr;
int_method_label = ['Using ' int_method ...
    ' interpolation before filtering in the time domain'];

switch int_method
    case 'none'
        disp(int_method_label)
    
    case 'linear'
        % calculate linear interpolation on high motion data using 
        % low-motion bold data
        temp=interp1(ix_in,Rr(ix_in,:),ix_ol);
        Rr_int(ix_ol,:)=temp;
        mean_Rr=mean(Rr(ix_in,:),1);
        
        if ~th(end) % if the last point has movement, use the mean
            Rr_int(end,:)=mean_Rr;
        end
        
        find_nans=mean(temp,2);
        find_nans=find(isnan(find_nans));
        Rr_int(ix_ol(find_nans),:)=repmat(mean_Rr,numel(find_nans),1);
        
        disp(int_method_label)
    otherwise
        disp('Method not found, using none')
end

%% band pass filter
F = 1/TR;
Ny = F/2;
BW_Hz = [lp_Hz hp_Hz];
BW_N = BW_Hz/Ny;
[b, a] = butter(ceil(bp_order/2),BW_N);
Y = filtfilt(b,a,Rr_int);
DVAR_post_all = dvars_from_cifti(Y);

% Quick visualization of results
%% remake cifti
dlmwrite(cifti_txt_path, Y','delimiter' ,' ');

cmd = [path_wb_c ' -cifti-convert -from-text ' cifti_txt_path ' ' path_cii ...
       ' ' result_dir filesep FNL_preproc_CIFTI_basename '.dtseries.nii'];

disp(cmd)
system(cmd);
system(['rm -f ' cifti_txt_path]);

%% Create grayplot
% adding nans to identify the begining of the sequence
nan_FD = [nan; FD];
nan_DVAR_pre_reg = [nan' DVAR_pre_reg];
nan_DVAR_post_reg = [nan' DVAR_post_reg];
nan_DVAR_post_all = [nan' DVAR_post_all];
%save('/mnt/max/home/andersp/workspace2.mat'

%load('/mnt/max/home/andersp/workspace1.mat')
fig_fMRI_QA(fMRIName,nan_FD, Xd, nan_DVAR_pre_reg, nan_DVAR_post_reg, ...
            nan_DVAR_post_all,path_ex_sum);
postreg_fig_fMRI_QA(fMRIName,nan_FD, Rr, nan_DVAR_pre_reg, ...
    nan_DVAR_post_reg, nan_DVAR_post_all,path_ex_sum)

% save data to be concatenated
temp_mat = [path_ex_sum filesep 'temp_grayplotdata_' fMRIName '.mat'];
save(temp_mat,'nan_FD');
save(temp_mat,'nan_DVAR_pre_reg','-append');
save(temp_mat,'nan_DVAR_post_reg','-append');
save(temp_mat,'nan_DVAR_post_all','-append');
save(temp_mat,'Xd','-append');
save(temp_mat,'Rr','-append');

%% Make motion_numbers.txt
filename = [result_dir filesep motion_filename];

is = 16;
formatSpec = ['%' num2str(is) 's %' num2str(is) 's %' num2str(is) 's %' ...
    num2str(is) 's %' num2str(is) 's %' num2str(is) 's %' num2str(is) ...
    's %' num2str(is) 's %' num2str(is) 's %' num2str(is) 's %' ...
    num2str(is) 's\n'];
fileID = fopen(filename,'w');
fprintf(fileID, formatSpec,'#difference', 'diff_x', 'diff_y', ...
    'diff_z', 'diff_arc_len_x', 'diff_arc_len_y', ...
    'diff_arc_len_z','FD','DVAR_pre_reg','DVAR_post_reg', 'DVAR_post_all');

ix{min(r,c)-1} = [];
dMR = diff(MR);
num_fmt = '%16.6f';
formatSpec = ['%' num2str(is) 's ' num_fmt ' ' num_fmt ' ' num_fmt ' ' ...
    num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt ' ' num_fmt ...
    ' ' num_fmt '\n'];
for i=2:min(r,c)
    ix{i-1} = [' ', num2str(i,'%04d') ' - ' num2str(i-1,'%04d') ' '];
    fprintf(fileID,formatSpec, ix{i-1}, dMR(i-1,:), FD(i-1), ...
        DVAR_pre_reg(i-1), DVAR_post_reg(i-1), DVAR_post_all(i-1));
end

fclose(fileID);
exit
