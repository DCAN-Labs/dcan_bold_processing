function [subject_survival_list, subject_survival, missing_data] = power_motion_survival_plot(group_folder_path,subject_list_file)
%Displays a colour map of the number of seconds remaining at various FD thresholds 

% Inputs:
% group_folder_path     Full filepath to where your timecourses are located
% subject_list_file     The name of your subject list .txt file

% Outputs:
% subject_list          A cell of subject names
% subject_survival      A matrix of time remaining at a given FD threshold
%                       for each subject (dimenions n by 51 where 51 is each FD threshold from 0
%                       to 0.5 in 0.1 increments) Subject order is the same
%                       as in subject_list.

% read in subject list
fid=fopen(subject_list_file);
vc_names=textscan(fid,'%s');
fclose(fid);

subject_list = vc_names{1};
subject_count = length(subject_list);

subject_survival = [];
missing_data = {};
subject_survival_list={};
row_count=1;

% loop through subjects
for i = 1:subject_count % for each subject
    sub_name = subject_list{i};
    disp(['Starting analyzing visit ' sub_name])
    subject_dir = [group_folder_path filesep sub_name];
    
    FD_path = [subject_dir filesep 'motion' filesep 'power_2014_motion.mat'];
    
    try
        load(FD_path);
    % for every subject get remaining seconds at each FD threshold
        for j = 1:size(motion_data,2);
            subject_survival(row_count,j) = motion_data{j,21}.remaining_seconds/60;
        end
        subject_survival_list=[subject_survival_list;{sub_name}];
        row_count=row_count+1;
        
    catch FD_mat_exception
        disp(['Excluding ' sub_name ' (Missing or broken  power_2014_motion.mat link): ' FD_path])
         missing_data = [missing_data; {sub_name}];
    end
    

end

% plot result
figure1 = figure('Colormap',...
    [1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;0.6875 1 0.3125;0.589285731315613 1 0.41071429848671;0.491071432828903 1 0.508928596973419;0.392857134342194 1 0.607142865657806;0.294642865657806 1 0.705357134342194;0.196428567171097 1 0.803571403026581;0.0982142835855484 1 0.901785731315613;0 1 1;0 0.9375 1;0 0.875 1;0 0.8125 1;0 0.75 1;0 0.6875 1;0 0.625 1;0 0.5625 1;0 0.5 1;0 0.4375 1;0 0.375 1;0 0.3125 1;0 0.25 1;0 0.1875 1;0 0.125 1;0 0.0625 1;0 0 1;0 0 0.9375;0 0 0.875;0 0 0.8125;0 0 0.75;0 0 0.6875;0 0 0.625;0 0 0.5625]);
imagesc(subject_survival)
set(gcf,'Position',[677 490 560 420])

xlabel('FD Threshold')
ylabel('Subject Number')
set(gca,'XTick',[6,11,16,21,26,31,36,41,46,51],'XTickLabel',{'0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
colorbar
