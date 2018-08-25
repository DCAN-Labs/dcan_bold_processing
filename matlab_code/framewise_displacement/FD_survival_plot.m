function [subject_survival_list, subject_survival, missing_data] = FD_survival_plot(group_folder_path,subject_list_file)
%Displays a colour map of the number of seconds remaining at various FD thresholds

% Inputs:
% group_folder_path     Full filepath to where your subject visit folders
%                       are located with timecourses and FD mat files.
% subject_list_file     The full path to your subject list .txt file
% 
% Outputs:
% subject_survival_list A cell of subject names who have valid FD data  (Changed by LKV 7/16/2014)
% subject_survival      A matrix of time remaining at a given FD threshold
%                       for each subject (dimenions n by 51 where 51 is each FD threshold from 0
%                       to 0.5 in 0.1 increments) Subject order is the same
%                       as in subject_list.  (Added by LKV 7/16/2014)
% missing_data          A cell of subject names whose FD data is missing
%                       (Added by LKV 7/16/2014)

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
    
    FD_path = [subject_dir filesep 'motion' filesep 'FD.mat'];
    
    try
        load(FD_path);
        % for every subject get remaining seconds at each FD threshold
        for j = 1:size(FD_data,2)
            subject_survival(row_count,j) = FD_data{j}.remaining_seconds/60;
        end
        subject_survival_list=[subject_survival_list;{sub_name}];
        row_count=row_count+1;
        
    catch FD_mat_exception
        disp(['Excluding ' sub_name ' (Missing or broken FD.mat link): ' FD_path])
        missing_data = [missing_data; {sub_name}];
    end
    
    
    
end

% plot result
figure;
imagesc(subject_survival)
set(gcf,'Position',[677 490 560 420])
mappers = colormap;
reverse_order = length(mappers):-1:1;
colormap(mappers(reverse_order',:))

xlabel('FD Threshold')
ylabel('Subject Number')
set(gca,'XTick',[6,11,16,21,26,31,36,41,46,51],'XTickLabel',{'0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
colorbar
