%% Figure4: SNR Models
% Figure layout
clear

% Set dots per inch
dpi = 96;

figure('units','pixels','Position',[2*dpi 0 8.5*dpi 11*dpi])

% From bottom to top
V = [1.5 3.5 5.5 7.5];

% From left to right
H = [1.5 4.1];

% Layout of axes
pos      = [H(1)  V(4)  1.5 1.5;
            H(2)  V(4)  1.5 1.5;
            
            H(1)  V(3)  1.5 1.5;
            H(2)  V(3)  1.5 1.5;
            
            H(1)  V(2)  1.5 1.5;
            H(2)  V(2)  1.5 1.5;
            
            H(1)  V(1)  1.5 1.5;
            H(2)  V(1)  1.5 1.5].*dpi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR Models: Increasing Signal and Decreasing Noise Axes 1 - 4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the relevant files
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
file1 = 'Figure5_ICMS_PSTH_SNR_bin10_win500.mat';
file2 = 'Learning_ICMS_SNR_models.mat';
load(file1)
load(file2)

% Axes 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of Unit excited by ICMS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line = -500:10:490;
axes('units','pixels','position',pos(1,:)),bar(time_line,R29_ICMS_psth_2(:,12),'k','barwidth',1)
    xlim([-500 500])
    ylabel('spikes/s','fontname','arial','fontsize',10)
    title([{'SNR Models:'}; {'ICMS Signal'}],'fontname','arial','fontsize',10,'fontweight','bold')

% Axes 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of SNR signal Increase Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('units','pixels','position',pos(2,:)),bar(Smeans,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       line([1:3; 1:3],[Smeans - Sstds; Smeans + Sstds],'color','k','linewidth',3)
       set(gca,'xticklabel','1|2|3')
       ylabel('z-score')
       title([{'Group Average by'}; {'Behavioral Stage'}],'fontname','arial','fontsize',10,'fontweight','bold')
       set(gca,'fontname','arial','fontsize',10)
       
% Axes 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of Unit inhibited by ICMS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line = -500:10:490;
axes('units','pixels','position',pos(3,:)),bar(time_line,V7_ICMS_psth_1(:,11),'k','barwidth',1)
    xlim([-500 500])
    ylabel('spikes/s','fontname','arial','fontsize',10)
    title('Internal Noise','fontname','arial','fontsize',10,'fontweight','bold')
    
% Axes 4
axes('units','pixels','position',pos(4,:)),bar(NImeans,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       line([1:3; 1:3],[NImeans - NIstds; NImeans + NIstds],'color','k','linewidth',3)     
       set(gca,'xticklabel','1|2|3')
       ylabel('t-value')
       set(gca,'fontname','arial','fontsize',10)
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR Models: Increasing Whisker System Sensitivity Axes 5 - 8  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the relevant files
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
file = 'Figure5_whisk_mod_depth_test_bin20_win1000.mat';
load(file)

% Axes 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of Unit excited by Whisker Stimulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line = -1000:20:990;
axes('units','pixels','position',pos(5,:)),bar(time_line,R29_whisk_psth_3(:,2),'k','barwidth',1)
    xlim([-1000 1000])
    ylabel('spikes/s','fontname','arial','fontsize',10)
    title('External Drive: Excitatory','fontname','arial','fontsize',10,'fontweight','bold')

% Axes 6
axes('units','pixels','position',pos(6,:)),bar(NEmeans_e,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       line([1:3; 1:3],[NEmeans_e - NEstds_e; NEmeans_e + NEstds_e],'color','k','linewidth',3)
       set(gca,'xticklabel','1|2|3')
       ylabel('z-score')
       set(gca,'fontname','arial','fontsize',10)

% Axes 7
axes('units','pixels','position',pos(7,:)),bar(time_line,V8_whisk_psth_2(:,7),'k','barwidth',1)
    xlim([-1000 1000])
    ylabel('spikes/s','fontname','arial','fontsize',10)
    xlabel('Time (msec)')
    title('External Drive: Inhibitory','fontname','arial','fontsize',10,'fontweight','bold')

% Axes 8
axes('units','pixels','position',pos(8,:)),bar(NEmeans_i,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       line([1:3; 1:3],[NEmeans_i - NEstds_i; NEmeans_i + NEstds_i],'color','k','linewidth',3)
       set(gca,'xticklabel','1|2|3')
       ylabel('z-score')
       set(gca,'fontname','arial','fontsize',10)

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')

set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 11])
