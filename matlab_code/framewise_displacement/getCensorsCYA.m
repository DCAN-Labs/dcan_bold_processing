function getCensorsCYA(working_dir,project_dir,fd_thresh,dvar_thresh)

%% Script to extract censor vectors for fMRI analyses [Created by Marc Rudolph. Oregon Health & Science University 2014]
% Update 9/30/2014 MR: Now can select which motion correction .mat file to use (Current FD.mat (FD Only) or Power2014.mat (FD + DVARS)

% Example function calls: 
% FD only:     getCensorsCYA('/group_shares/NAGEL_LAB/staff/gaby/flanker/CYA','/group_shares/NAGEL_LAB/CYA/analyses_v2/Flanker',.3)
% FD + DVARs:  getCensorsCYA('/group_shares/NAGEL_LAB/staff/gaby/flanker/CYA','/group_shares/NAGEL_LAB/CYA/analyses_v2/Flanker',.3,20)

% addpath(genpath('/group_shares/FAIR_LAB/usr/Matlab_Scripts'));
% addpath(genpath('/group_shares/PSYCH/code/development/utilities/framewise_displacement'));


%% Setup directories
currdir = pwd;
cd(working_dir);

%% Setup Subjects
subject_info_file = uigetfile('*.txt', 'Select a subject list TEXT (.txt) file');
fid=fopen(subject_info_file);
vc_names=textscan(fid,'%s');
fclose(fid);
subs=length(vc_names{1,1});

%%Setup thresholds
fd_thresh = (fd_thresh*100)+1;
if exist('dvar_thresh','var') == 1
dvar_thresh = dvar_thresh+1;
end

%% Create Sub Censor Files
for ki=1:subs
    
    % Setup subjects
    sub_name=vc_names{1,1}{ki,1};
    subj = strsplit(sub_name,'+');      % Split string to spearate scan & ID (**matlab13a required**)
    sub = subj{1}; scan = subj{2};      % Subject ID & Scandate
    
    % Set filepaths
    cd(project_dir);                    % Where are your processed files located - lets go there
    pipeline = dir('F*');               % Load the motion .mat file
    matdir = [pwd '/' pipeline.name '/' sub_name '/motion'];
    
    
    %Setup censors for FD or FD + DVARs
    censors = {};
        if exist('dvar_thresh','var') == 1
            display('Using FD + DVARs')
            matfile = [matdir '/power_2014_motion.mat']; load(matfile);
            censored_frames = motion_data{fd_thresh,dvar_thresh}.frame_removal;
        else
            display('Using FD ONLY')
            matfile = [matdir '/FD.mat']; load(matfile);
            censored_frames = FD_data{1,fd_thresh}.frame_removal;
        end
    censors = censored_frames;
   
    
    % Save out vector in the current subjects folder
    cd(working_dir);
    
        if exist(sub,'dir') == 0
            mkdir(sub);
            cd(sub);
        else
            disp('subject directory exists - just writing file');
            cd(sub);
        end
        
    dlmwrite('censors.txt',censors);

    % Return the working directory
    cd(working_dir);
   
    
end