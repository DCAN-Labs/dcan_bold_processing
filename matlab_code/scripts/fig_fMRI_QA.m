function fig_fMRI_QA(fMRIName,nan_FD, Xd, nan_DVAR_pre_reg, nan_DVAR_post_reg, nan_DVAR_post_all,path_ex_sum)
%% This function shows gray plots and some extra info required for quick visual inspection of bold data
% Oscar Miranda-Dominguez

%silence warnings
warning('off', 'all')


[r c]=size(Xd);
%% settings

FDvlow_th = 0.05;
FDlow_th  = 0.10;
FDmid_th  = 0.20;
FDhigh_th = 0.50;
FDmax_th  = 0.80;
FD_ytick_res = 0.1;
WBmin = -400;
WBmax = 800;
DVARmax = 500;
DVAR_ytick_res = 100;

n_bins=21;

my_color=[
    228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0;
    0,127,0;
    150,150,150;
    220,220,0
    ]/255;

red_ish       = my_color(1,:);
blue_ish      = my_color(2,:);
litegreen_ish = my_color(3,:);
purple_ish    = my_color(4,:);
orange_ish    = my_color(5,:);
darkgreen_ish = my_color(6,:);
gray_ish      = my_color(7,:);
yellow_ish    = my_color(8,:);

lw=1; % plot line widths
ms=7; % marker sizes
%%


fig_wide=18;

lm=1.5;
pre_h=0;
bw_h=.2;
pw=8;% panel width
pos_h=4/fig_wide;
rm=1.5;
sum_x=lm+pre_h+bw_h+pw+bw_h+pos_h+rm;

% Normalize x
lm=lm/sum_x;
pre_h=pre_h/sum_x;
bw_h=bw_h/sum_x;
pw=pw/sum_x;
pos_h=pos_h/sum_x;
rm=rm/sum_x;

% calculate the required hight in cm

bm=1.5;
bw_v=.15;
ph=3;
tm=.75;

sum_y=bm+(ph+bw_v)+(2*ph+bw_v)+(ph+bw_v)+(ph)+tm;
bm=bm/sum_y;
bw_v=bw_v/sum_y;
ph=ph/sum_y;
tm=tm/sum_y;
fig_tall=sum_y;

%%
h = figure('Visible','on',...
    'units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[fig_tall fig_wide],...
    'name',fMRIName,...
    'Position',[8 1 fig_wide fig_tall],...
    'PaperPosition',[8 1 fig_wide fig_tall],...
    'color',[1 1 1]);
set(h, 'PaperUnits','centimeters')
set(h,'PaperSize',[fig_tall fig_wide] )
set(h,'PaperPosition',[8 1 fig_wide fig_tall] )


fd_panel=subplot('position',[(lm+pre_h+bw_h) bm pw ph]);
gr_panel=subplot('position',[(lm+pre_h+bw_h) (bm+ph+bw_v) pw 2*ph]);
gs_panel=subplot('position',[(lm+pre_h+bw_h) (bm+3*ph+2*bw_v) pw ph]);
dv_panel=subplot('position',[(lm+pre_h+bw_h) (bm+4*ph+3*bw_v) pw ph]);

pos_fd_panel=subplot('position',[(lm+pre_h+2*bw_h+pw) bm pos_h ph]);
pos_gr_panel=subplot('position',[(lm+pre_h+2*bw_h+pw) (bm+ph+bw_v) pos_h 2*ph]);
pos_gs_panel=subplot('position',[(lm+pre_h+2*bw_h+pw) (bm+3*ph+2*bw_v) pos_h ph]);
pos_dv_panel=subplot('position',[(lm+pre_h+2*bw_h+pw) (bm+4*ph+3*bw_v) pos_h ph]);


%%
subplot(fd_panel)

%replace all nan with 0
%nan_FD(isnan(nan_FD)) = 0;

undervlowFD_ix=find(nan_FD<FDvlow_th);
vlowFD_ix=find(nan_FD>=FDvlow_th & nan_FD<FDlow_th);
lowFD_ix=find(nan_FD>=FDlow_th & nan_FD<FDmid_th);
midFD_ix=find(nan_FD>=FDmid_th & nan_FD<FDhigh_th);
highFD_ix=find(nan_FD>=FDhigh_th);
overmaxFD_ix=find(nan_FD>FDmax_th);
nan_ix=find(isnan(nan_FD));


plot(1:r,nan_FD,'k','linewidth',lw)
hold all
%add threshold lines at 0.2 and 0.5 and mark the points - AP 20170303
plot([1 r],FDvlow_th*[1 1],'-', 'color', gray_ish,'linewidth',lw)
plot([1 r],FDlow_th*[1 1],'-', 'color', darkgreen_ish,'linewidth',lw)
plot([1 r],FDmid_th*[1 1],'-', 'color', orange_ish,'linewidth',lw)
plot([1 r],FDhigh_th*[1 1],'r-','linewidth',lw)
plot(vlowFD_ix,nan_FD(vlowFD_ix),'.','MarkerSize',ms,'color', gray_ish);
plot(lowFD_ix,nan_FD(lowFD_ix),'.','MarkerSize',ms,'color', darkgreen_ish);
plot(midFD_ix,nan_FD(midFD_ix),'.','MarkerSize',ms,'color', orange_ish);
plot(highFD_ix,nan_FD(highFD_ix),'.r','MarkerSize',ms);
%markers for the top of the FD plot
% plot(underlowFD_ix,FDmax_th,'.','MarkerSize',ms,'color', darkgreen_ish);
if ~isempty(vlowFD_ix)
    plot(vlowFD_ix,FDmax_th,'.','MarkerSize',ms,'color', gray_ish);
end
if ~isempty(lowFD_ix)
    plot(lowFD_ix,FDmax_th,'.','MarkerSize',ms,'color', darkgreen_ish);
end
if ~isempty(midFD_ix)
    plot(midFD_ix,FDmax_th,'.','MarkerSize',ms,'color', orange_ish);
end
if ~isempty(highFD_ix)
    plot(highFD_ix,FDmax_th,'.r','MarkerSize',ms);
end
if ~isempty(overmaxFD_ix)
    plot(overmaxFD_ix,FDmax_th,'.r','MarkerSize',ms);
end

for i=2:length(nan_ix)
    line(nan_ix(i)*[1 1],[0 FDmax_th],'color',[1 1 1]*.7);
end

hold off
ylabel('FD')
xlabel('Frame')
axis tight
ylim([0 FDmax_th]) 
grid off
set(gca,'ytick',0:FD_ytick_res:FDmax_th);
set(gca,'yticklabel',num2str(get(gca,'ytick')','%4.1f'))
xl=xlim;

% grayplot modifications - AP 20171003
cl=[1 99]; % percentiles to be shown
% cl=[0 100]; % percentiles to be shown

%clg=prctile(Xd(:),cl);
clg=[-600 600];

subplot(gr_panel)
sampled_Xd=Xd(:,1:10:length(Xd));
%imagesc((1:r)-.5,1:c,sampled_Xd')
imagesc((1:r)-.5,1:c,Xd',clg)
colormap(gray)

hold on
plot([1 r],[c*(1-(91282-2*32492)/91282) c*(1-(91282-2*32492)/91282)],'--k','linewidth',lw)
hold off

set(gca,'ytick',[32492 2*32492])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
ylabel({num2str(c), ' grayordinates'})
xlim(xl)
text(xl(2)+0.1,c*32492/91282,'CORTEX',...
    'VerticalAlignment','baseline','fontsize',8,'color','k','Rotation',-90,'HorizontalAlignment','center')
text(xl(2)+0.1,c*(1-(91282-2*32492)/91282/2),'SUBCORT',...
    'VerticalAlignment','baseline','fontsize',8,'color','k','Rotation',-90,'HorizontalAlignment','center')
text(xl(2)+60,c*(1-(91282-2*32492)/91282/2)+11000,num2str(clg(1)),'VerticalAlignment','baseline','fontsize',8,'color','k','HorizontalAlignment','center')
text(xl(2)+60,c*32492/91282-28000,num2str(clg(2)),'VerticalAlignment','baseline','fontsize',8,'color','k','HorizontalAlignment','center')


subplot(gs_panel)
plot(mean(Xd,2),'color',red_ish,'linewidth',lw)
hold all
plot(std(Xd,[],2),'color',blue_ish,'linewidth',lw)
hold off
set(gca,'xticklabel',[])
%remove autoscaling of FD plot. Set y-axis range to -400,1000 - AP 20170303
%axis tight
axis([xl WBmin WBmax])
ylabel('WB')
grid off

hold all
for i=2:length(nan_ix)
    line(nan_ix(i)*[1 1],ylim,'color',[1 1 1]*.7)
end
hold off

subplot(dv_panel)
plot(1:r,nan_DVAR_pre_reg,'color',litegreen_ish,'linewidth',lw);
hold all
plot(1:r,nan_DVAR_post_reg,'color',purple_ish,'linewidth',lw);
plot(1:r,nan_DVAR_post_all,'color',orange_ish,'linewidth',lw);
hold off
ylabel('DVAR')
%remove autoscaling of DVARS plot. Set y-axis range - AP 20170303
axis([xl 0 DVARmax])
grid off
set(gca,'xticklabel',[])
set(gca,'ytick',0:DVAR_ytick_res:DVARmax)
xlim(xl)
hold all
for i=2:length(nan_ix)
    line(nan_ix(i)*[1 1],ylim,'color',[1 1 1]*.7)
end
hold off
title(fMRIName)

%%
subplot(pos_fd_panel)
x_hist=linspace(0,FDmax_th,n_bins);
y_hist=histc(nan_FD,x_hist);
barh(x_hist,y_hist)
hold all
%include threshold lines - AP 20170603
plot(xlim,FDvlow_th*[1 1],'-','color',gray_ish,'linewidth',lw)
plot(xlim,FDlow_th*[1 1],'-','color',darkgreen_ish,'linewidth',lw)
plot(xlim,FDmid_th*[1 1],'-','color',orange_ish,'linewidth',lw)
plot(xlim,FDhigh_th*[1 1],'r-','linewidth',lw)
hold off
ylim([0 FDmax_th]) 
set(gca,'ytick',0:FD_ytick_res:FDmax_th);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
set(gca, 'xtick', [])

xl=xlim;

bottom_arrow_loc = FDvlow_th/2;
text(xl(2),bottom_arrow_loc,'\leftarrow','fontsize',8,'color','k')
text(xl(2),bottom_arrow_loc,[ '    ' num2str(length(undervlowFD_ix)) ' frames'],...
    'VerticalAlignment','middle','fontsize',6,'color','k')

vlow_arrow_loc = (FDlow_th-FDvlow_th)/2 + FDvlow_th;
text(xl(2),vlow_arrow_loc,'\leftarrow','fontsize',8,'color',gray_ish)
text(xl(2),vlow_arrow_loc,[ '    ' num2str(length(undervlowFD_ix)+length(vlowFD_ix)) ' frames'],...
    'VerticalAlignment','middle','fontsize',6,'color',gray_ish)

low_arrow_loc = (FDmid_th-FDlow_th)/2 + FDlow_th;
text(xl(2),low_arrow_loc,'\leftarrow','fontsize',8,'color',darkgreen_ish)
text(xl(2),low_arrow_loc,[ '    ' num2str(length(undervlowFD_ix)+length(vlowFD_ix)+length(lowFD_ix)) ' frames'],...
    'VerticalAlignment','middle','fontsize',6,'color',darkgreen_ish)

mid_arrow_loc = (FDhigh_th-FDmid_th)/2 + FDmid_th;
text(xl(2),mid_arrow_loc,'\leftarrow','fontsize',8,'color',orange_ish)
text(xl(2),mid_arrow_loc,[ '    ' num2str(length(undervlowFD_ix)+length(vlowFD_ix)+length(lowFD_ix)+length(midFD_ix)) ' frames'],...
    'VerticalAlignment','middle','fontsize',6,'color',orange_ish)

high_arrow_loc = (FDmax_th-FDhigh_th)/2 + FDhigh_th;
text(xl(2),high_arrow_loc,'\leftarrow','fontsize',8,'color','r')
text(xl(2),high_arrow_loc,[ '    ' num2str(length(undervlowFD_ix)+length(vlowFD_ix)+length(lowFD_ix)+length(midFD_ix)+length(highFD_ix)) ' frames'],...
    'VerticalAlignment','middle',...
    'fontsize',6,'color','r')

%hack a scale using contourf, since imagesc errors on octave
subplot(pos_gr_panel)
xl=xlim;
Xd_min = min(sampled_Xd(:));
Xd_max = max(sampled_Xd(:));
Xd_range = abs(Xd_max - Xd_min);
Xd_scale = [[Xd_min:Xd_max]',[Xd_min:Xd_max]'];
%[cv,ch] = contourf(Xd_scale,100);

local_scale=[clg(1):clg(2)];
if numel(local_scale)>2000
    local_scale=linspace(clg(1),clg(2),2000);
end
imagesc(local_scale')
%set(ch,'LineColor','none');
colormap(gray);
set(gca,'YTickLabel',[]);
set(gca, 'ytick', [])
set(gca, 'xtick', [])
set(gca, 'xticklabel', [])


subplot(pos_gs_panel)
yl=ylim;
yp=linspace(yl(1),yl(2),4);
text(0, yp(2),'Mean','color',red_ish)
text(0, yp(3),'Std','color',blue_ish)

axis off
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
set(gca, 'ytick', [])
set(gca, 'xtick', [])

subplot(pos_dv_panel)
yl=ylim;
yp=linspace(yl(1),yl(2),5);
text(0, yp(4),'Pre reg','color',litegreen_ish)
text(0, yp(3),'Post reg','color',purple_ish)
text(0, yp(2),'Post all','color',orange_ish)

set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off
png_file=[path_ex_sum '/DVARS_and_FD_' fMRIName '.png'];
print(h,png_file, '-dpng','-r300')
disp(png_file)
disp('DONE')
end
% exit
