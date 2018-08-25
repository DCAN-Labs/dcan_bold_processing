 %function create_FD_report(group, basedir, need_sublist,need_demo,cutoff,fd_thresh)
%function create_FD_report(group, basedir, need_sublist,need_demo,cutoff,fd_thresh,view_plots)

function create_FD_report(group, fd_thresh, cutoff)

%% Instructions - Expand to view:

% This script creates movement reports (FD numbers, Censors, List of Movers) - Steps:

% Go to your studies processed directory (this is basedir) & create a sublist (edit after to sort & include who you want).
% cd /group_shares/FAIR_LAB2/CYA/processed/ADHD-HumanYouth-OHSU

% slither -t -f --max 6 --min 6 -m '(.+)FD.mat'  . 'echo \g<0>'
% Take off the -t after checking the paths are correct

% Go to (cd) your project/personal directory

% Call the function as below examples:
% create_FD_report('ADHD_TestFD', '/group_shares/FAIR_LAB2/CYA/processed/ADHD-HumanYouth-OHSU',0,0,5, [2,3])
% create_FD_report('ADHD_TestFD', '/group_shares/FAIR_LAB2/CYA/processed/ADHD-HumanYouth-OHSU', 'sub_test.txt',5, [2,2.5,3])
% create_FD_report('MacArt_Task', '/scratch/marc_temp/MacArthur_TaskBased/', 'subs.txt',3,[2,3,4,5])

%% Setup directories & subjects (Manually - Ignore)

%MacArthur Task
% group = 'MacArt_Task_173';
% basedir = '/scratch/marc_temp/MacArthur_TaskBased/';
% fileloc = '/FAIRPRE10_TR2pt5_RAD50pt0_SKIPSEC5pt0/FCON10/21nFDPAR_noopt/FD.mat';
% subject_info_file = 'Macarthur_98Subs_CYA_minusBadReg_20140715_old.txt';
% groupID = 'macart';
% cutoff = 180;
% fd_thresh = 3;

%ADHD - CYA
% group = 'ADHD_MarchApril';
% basedir = '/group_shares/FAIR_LAB2/CYA/processed/ADHD-HumanYouth-OHSU/';
% basedir = currdirr;
% subject_info_file  = 'slither_sublist_FD.txt';

%% Set current (working) directory & processed directory of choice
currdir = pwd;
% cutoff = 3.96;
% cutoff = 180; %seconds

screen = get(0,'ScreenSize');
set(0, 'DefaultFigurePosition', [screen/2]);


dir_options = {'/group_shares/FAIR_LAB2/CYA/processed/ADHD-HumanYouth-OHSU'...
    '/group_shares/FAIR_ASD/CYA/processed/ASD-HumanYouth-OHSU'...
    '/group_shares/FAIR_LAB2/CYA/processed/ADHD-ASD-HumanInfant-OHSU'...
    '/group_shares/FAIR_LAB2/CYA/processed/Infant-HumanInfant-UCIrvine'...
    '/group_shares/FAIR_LAB2/CYA/processed/EEG-HumanYouth-OHSU'...
    '/group_shares/FAIR_LAB2/CYA/processed/Music-HumanAdult-Auckland_OHSU'...
    '/scratch/marc_temp/processed/MacArthur-HumanYouth_TaskBased'...
    'Not listed, I will choose the directory myself'...
     '/group_shares/FAIR_LAB2/CYA/processed/Infant-HumanInfant-UCIrvine/ONE_YEAR'};


        %'/scratch/marc_temp/MacArthur_TaskBased'... 

basedir = menu('Select your studies processed directory', dir_options);


if basedir == 7
    groupID = 'macart';
else if basedir == 4 | basedir == 9
        groupID = 'irvine';
    else
        groupID = basedir;
    end
end


if basedir ~= 8
    display(['Study selected:' dir_options{basedir}]);
else
    basedir = uigetdir('/group_shares');
    display(['Study selected:' basedir]);
end
basedir = dir_options{basedir}


demo_options = menu('Need DEMO information?','YES (ADHD Only)','NO (Default)');
if demo_options == 1
    need_demo = 1;
else
    need_demo = 0;
end

plot_options = menu('Plot Results?','YES','NO');
if plot_options == 1
    view_plots = 1;
else
    view_plots = 0;
end

sublist_options = menu('Need a subject list?', 'YES','NO');
if sublist_options == 1
    need_sublist = 1;
else
    need_sublist = 0;
end


%% Slither created subject list - Expand to view:

if need_sublist == 1
    
    cd (basedir)
    
    syscom = [' slither -f --max 6 --min 6 -m ''(.+)FD.mat'' . ''echo \g<0>'' >> mysubs_temp.txt ' ];
    command = [syscom]; system(command);
    
    movefile('mysubs_temp.txt', currdir);

    
    cd (currdir);
    
    find_and_replace('mysubs_temp.txt','\./',''); %Backslash needed for '.' in regular expressions
    newfile = [group '_FDreport_subjectlist.txt']
    movefile('mysubs_temp.txt',newfile);
    syscom = [' gedit ' newfile ' & ']; command = [syscom]; system(command);
    pause   %Click any button to continue
    disp('Paused for editing. Press any key to continue after closing your subject list');
    subject_info_file = newfile;
    
else
    
    disp('Subject defined user list');
    subject_info_file = uigetfile('*.txt', 'Select a subject list TEXT (.txt) file');
    
end


fid=fopen(subject_info_file);
vc_names=textscan(fid,'%s');
fclose(fid);
subs=length(vc_names{1,1});

%% -or- Use All Subs in dir - Expand to view:
% subs=dir('1*');
% subs = {subs.name}';

%% If caselist info needed! - Expand to view:

if need_demo == 1
     disp('Loading Case Information');
    create_caselist_report(group, basedir, subject_info_file);
    case_dir = dir(['Case_info*' '.mat']);
    load(case_dir.name);
    cell_info = [];
else
    disp('no caselist information required');
end

%% Global setting for viewing plots or not

if view_plots == 1
    disp('show plots');
    set(gcf,'Visible','on')
    set(0,'DefaultFigureVisible','on');
    n = 5; 
    x = rand(n);
    col = jet(6); % or some other colormap
else
    set(gcf,'Visible','off')              % turns current figure "off"
    set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
end




%% FD report for multiple thresholds : Initialize all variables of interest
fd_thresh = [fd_thresh]; %fd_thresh = [2,3];  % 2 = .2FD; 3 = .3FD

subjects = {}; subjectIDs = {}; scandate = {}; subject_survival = []; all_subs=[];

FD_Report = []; FD_Sec_Remaining = []; FD_Min_Remaining = []; Frame_Censor_Report = [];  Frame_Censors = [];
    Movers={}; Movers_List = []; Movers_Report = []; 
Full_report =[]; Full_Report_Matrix = []; Full_Report_Matrix2 =[];
Labels = {}; Final_Labels={}; Final_Labels2 = [];
scanonly = {};  finalsubs = []; finalscans = [];

%% Begin
for mv=1:length(fd_thresh)

    thresh = ((fd_thresh(mv)*10)+1);    % CYA: Load FD Thresh from correct cell in matrix;

    for ki=1:subs
        
        sub_name=vc_names{1,1}{ki,1};       %Split string - matlab13a required***
        subj = strsplit(sub_name,'/');
        sub = subj{1};
        scan = subj{2};
        scanonly = strsplit(scan,'-');
        scandate = scanonly{1};
        subjects{end+1} = sub;
 

        %Merge ID's
        ID_temp = char(subj(1));
        
        if groupID == 'macart'
            ID = strrep(ID_temp,'SLA_','');
            ID = strrep(ID,'SC0_','');
            ID = strrep(ID,'SC0','');
            ID = strrep(ID,'_','');
        else if groupID == 'irvine'
                ID = strrep(ID_temp,'BUSS_','');
            else
                ID = strrep(ID_temp,'-','');
                ID = strrep(ID,'W1',''); ID = strrep(ID,'w2',''); ID = strrep(ID,'w3','');
            end
        end


        subjectIDs{end+1} = ID;
        finalsubs(ki,1) = str2double(ID);
        finalscans(ki,1) = str2double(scandate);
 
        matfile = sub_name;
        cd (basedir)
        load(matfile);
        
%         ki = 1; mv = 1;
%         thresh = [3];
        FD_Sec_Remaining(ki,mv) = FD_data{1,thresh}.remaining_seconds;
        FD_Min_Remaining(ki,mv) = ((FD_data{1,thresh}.remaining_seconds)/60);
        Frames_Remaining(ki,mv) = FD_data{1,thresh}.remaining_frame_count;
        Percent_Frames_Remaining(ki,mv) = ((Frames_Remaining(ki,mv) / FD_data{1,1}.total_frame_count)*100);
        
        
        %% Get Censors
        
        %getCensorsCYAPP_mvmtreport
        
        % Frame_length = length(FD_data{1,thresh}.frame_removal);
        % Frame_Censors = zeros(Frame_length,1);
        % Frame_Censors = FD_data{1,thresh}.frame_removal;
        % Frame_Censor_Report(:,ki) = Frame_Censors;
        
        % censors = {};
        % censored_frames = FD_data{1,mv}.frame_removal;
        % censors = [censors censored_frames];
        % dlmwrite('censors.txt',censors);
        
        %% Make column vectors
        %         FD_Min_Remaining = FD_Min_Remaining.';
        %         %FD_vector(size(FD_Report),ki) = FD_Report.';
        %         FramesRemaining = Frames_Remaining.';
        %         Percent_Remaining = Percent_Frames_Remaining.';
        
        %%  Movers Report
        %cutoff = (((FD_data{1,1}.total_frame_count * FD_data{1,thresh}.epi_TR)/60)*.5);
      
       if FD_Sec_Remaining(ki,mv) < cutoff
        Movers_List(ki,1) = finalsubs(ki,1);
        Movers_List(ki,2) = finalscans(ki,1);
        end
        

  %% Populate case info
        if need_demo == 1
            for c=1:6
                cell_info(ki,c) = Case_info{2,1}(ki,c);
            end
        else
        end
        
        %% Incorporating Survival Plot (Will use function call later....)
        for j = 1:size(FD_data,2)
            if FD_data{j}.remaining_seconds/60 ~= 0
                subject_survival(ki,j) = FD_data{j}.remaining_seconds/60;
            else
                subject_survival(ki,j) = 0;
            end
        end
        
        
end   % End subjects loop
    
    
    %% CD back to working directory & save
    
    cd (currdir);
    
    Full_Report = {FD_Sec_Remaining(:,mv), FD_Min_Remaining(:,mv), Frames_Remaining(:,mv), Percent_Frames_Remaining(:,mv) };
    Full_Report_Matrix(:,:,mv)  = [FD_Sec_Remaining(:,mv) FD_Min_Remaining(:,mv) Frames_Remaining(:,mv) Percent_Frames_Remaining(:,mv)];
    %CV = [ {['FD_' num2str(fd_thresh(mv)) '_Seconds_Remaining'],['FD_' num2str(fd_thresh(mv)) '_Minutes_Remaining'],['FD_' num2str(fd_thresh(mv)) '_Frames_Remaining'] ,['FD_' num2str(fd_thresh(mv)) '_Percent_Frames_Remaining']};{Full_Report{1}, Full_Report{2}, Full_Report{3}, Full_Report{4}}];
    
    % Generate Labels for each threshold
    Labels{mv} = {['FD pt' num2str(fd_thresh(mv)) 'mm Seconds Remaining'], ['FD pt' num2str(fd_thresh(mv)) 'mm Minutes Remaining'], ['FD pt' num2str(fd_thresh(mv)) 'mm Total Frames Remaining'], ['FD pt' num2str(fd_thresh(mv)) 'mm Percent Frames Remaining']};
    Final_Labels(1,:,mv) = Labels{mv};
    
    
%% Output list of movers for each thresh ?
% Movers_Report(:,:) = [Movers_List];
% Movers_Report2 = [];
% m = length(Movers_List);
% Movers_Report2 = [reshape(Movers_Report,m,[])]; 
% 
% indices = find(Movers_Report2(:,2)==0);
% Movers_Report2(indices,:) = []; 
% total_movers = length(Movers_Report2);
% display([num2str(total_movers) ' subjects with less than ' num2str(cutoff) ' seconds remaining at pt ' num2str(fd_thresh(mv)) ' FD']);
% 
%   
%     fid = fopen([group '_movers_list_cuttoff_seconds_' num2str(cutoff) '_at_pt' num2str(fd_thresh(mv)) 'FD.txt'],'w') ;     % Open new/blank txt file
%     fprintf(fid,['The Movers \n']);
%     fclose(fid);
%     dlmwrite([group '_movers_list_cuttoff_seconds_' num2str(cutoff) '_at_pt' num2str(fd_thresh(mv)) 'FD.txt'],Movers_Report2, '-append', 'delimiter', '\t','precision', 8);   % Write data to text file - below headers
%     disp('List of movers generated at specified cutoff');
%     
% clear m indices total_movers  Movers_Report Movers_List Movers_Report2
% Movers_Report = []; Movers_List = []; Movers_Report2 = [];
%     
end   % End threshold loop


%% Build final labeled matrices

    
Final_Labels = {reshape(Final_Labels,1,[])};
Full_Report_Matrix2 = [Full_Report_Matrix2 Full_Report_Matrix(:,:,:)];      %Full_Report_Matrix2 = [Full_Report_Matrix(:,:,1) Full_Report_Matrix(:,:,2)];
Full_Report_Matrix3 = reshape(Full_Report_Matrix2,subs,[]);

%Full_Report_Matrix3 = [Full_Report_Matrix2];

%% Output report with Caselist DEMO information
if need_demo == 1
    for f=1:4
        cell_info(ki,f) = Case_info{2,1}(ki,f);
    end
    
    cell_info = [cell_info Full_Report_Matrix3];
    
end

%% Manually Label
%Labels={'Mergeid','W1_C_Dteamdcsn_Dtadhd','W1_C_Dteamdcsn_Overallsb','Access_Adhd_Status','Child_Dob','Sex'...
%         'FD .2mm Sec Remaining', 'FD .2mm Min Remaining', 'FD .2mm Frames Remaining', 'FD .2mm Percent Frames Remaining'...
%         'FD .3mm Sec Remaining', 'FD .3mm Min Remaining', 'FD .3mm Frames Remaining', 'FD .3mm Percent Frames Remaining' };
%
%     Final_Labels =[];
%     for l=1:8
%         Final_Labels = [Final_Labels Labels{l} '\t'];
%     end

% Create output
%     fid = fopen('mvmnt_report.csv','wt') ;     % Open new/blank txt file
%     fprintf(fid,[Final_Labels '\n']);
%     fclose(fid);
%     dlmwrite('mvmnt_report.csv',[cell_info], '-append', 'delimiter', '\t','precision', 6);   % Write data to text file - below headers


%% Generate Final Labels

for l=1:length(Final_Labels{1,:})
    Final_Labels2 = [Final_Labels2 Final_Labels{1,1}{1,l} '\t'];
end

Final_Labels3 = ['Subject ID\tScan Date\t' Final_Labels2];

%% Create report (with or without demo info)

if need_demo == 1
    fid = fopen([group '_Mvmnt_Report_withdemo_' date '.csv'],'wt') ;     % Open new/blank txt file
    fprintf(fid,['Scan Date \t Mergeid \t W1_C_Dteamdcsn_Dtadhd \t W1_C_Dteamdcsn_Overallsb \t Access_Adhd_Status \t Child_Dob \t Sex \t' Final_Labels2 '\n']);
    fclose(fid);
    dlmwrite([group '_Mvmnt_Report_withdemo_' date '.csv'],[finalscans cell_info], '-append', 'delimiter', '\t','precision', 8);   % Write data to text file - below headers
    disp('Movement Report Generated');
    %delete(case_dir.name);
else
    fid = fopen([group '_Mvmnt_Report_nodemo_' date '.csv'],'wt') ;     % Open new/blank txt file
    fprintf(fid,[Final_Labels3 '\n']);
    fclose(fid);
    dlmwrite([group '_Mvmnt_Report_nodemo_' date '.csv'],[finalsubs finalscans Full_Report_Matrix3], '-append', 'delimiter', '\t','precision', 8);   % Write data to text file - below headers
    disp('Movement Report Generated ');
end

save([group '_Mvmnt_Report_' date '.mat'],'finalsubs', 'finalscans', 'Full_Report_Matrix3');

%%

%% Plot result (exta option)
if view_plots == 1

bar(Full_Report_Matrix3,'stacked','DisplayName','Full_Report_Matrix3')
title('FD Data : Min to Max Threshold')
xlims = get(gca,'XLim');
hold on
plot(xlims, [285 285], '-r')
plot(xlims, [570 570], '-g')
plot(xlims, [855 855], '-b')
plot(xlims, [1140 1140], '-r')
set(gcf,'Position',[677 490 560 420]);

figure
for i=1:length(fd_thresh)
    hold on
plot(Full_Report_Matrix3(:,[4*i]),'color',col((i),:));
end
xlims = get(gca,'XLim');
hold on
plot(xlims, [30 30], '-r')
plot(xlims, [50 50], '-g')
title('Percent Frames remaining')
set(gca,'XMinorTick','on');
set(gcf,'Position',[600 490 800 420]);
hold off


figure
for i=1:length(fd_thresh)
    hold on
    if i == 1
plot(Full_Report_Matrix3(:,i),'color',col((i),:));
    else
        plot(Full_Report_Matrix3(:,((i-1)*5)),'color',col((i),:));
    end
end
xlims = get(gca,'XLim');
hold on
plot(xlims, [180 180], '-r')
title('Seconds remaining')
set(gca,'XMinorTick','on');
set(gcf,'Position',[600 490 800 420]);
hold off


figure
imagesc(subject_survival);
set(gcf,'Position',[677 490 560 420]);
mappers = colormap;
reverse_order = length(mappers):-1:1;
colormap(mappers(reverse_order',:))

xlabel('FD Threshold');
ylabel('Subject Number');
title('FD Frames Remaining : All Theshold : Per Subject')
set(gca,'XTick',[6,11,16,21,26,31,36,41,46,51],'XTickLabel',{'0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'});
colorbar
   
end


%% The End
disp('Voila - Your very own movement report, with a pretty picture, all in one simple step!');
clear all

end