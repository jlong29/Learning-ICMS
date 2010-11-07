%% Figure: Tracking SNR models
clear
close all

% Set dots per inch
dpi = 96;

figure('units','pixels','Position',[0 0 8.5 11].*dpi)

% From left to right
H = [0.8 2.3 3.8 5.3 6.8];

% From bottom to top
V = [1.95 3.6 4.85];

% Layout of axes
pos      = [H(1)  V(3)  1 1;
            H(2)  V(3)  1 1;
            H(3)  V(3)  1 1;
            H(4)  V(3)  1 1;
            H(5)  V(3)  1 1;
            
            H(1)  V(2)  1 1;
            H(2)  V(2)  1 1;
            H(3)  V(2)  1 1;
            H(4)  V(2)  1 1;
            H(5)  V(2)  1 1
            
            H(1)  V(1)  1 1;
            H(2)  V(1)  1 1;
            H(3)  V(1)  1 1;
            H(4)  V(1)  1 1;
            H(5)  V(1)  1 1].*dpi;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR Models and Subject Behavior: Axes 1 - 10 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
% ensemble SNR model data
load('Figure5_Tracking_SNR_models.mat')
file_list = {R29s,V1s,V4s,V7s,V8s};

% ICMS: excitation: signal
means = {R29_means; V1_means; V4_means; V7_means; V8_means};
ses   = {R29_SEs;   V1_SEs;   V4_SEs;   V7_SEs;   V8_SEs};

% ICMS: inhibition: internal noise
meansC = {R29_meanb; V1_meanb; V4_meanb; V7_meanb; V8_meanb};
sesC   = {R29_SEb;   V1_SEb;   V4_SEb;   V7_SEb;   V8_SEb};

% Subjects' behavioral data
load('Subjects_Behavior.mat')

% Axes pointers
axs1   = 1:5;
axs2   = 6:10;
axs3   = 11:15;

for p = 1:size(pos,1)/3
    beh_sum = subject_Behavior{p};
    
    if p == 1
        metric = eval(['length(file_list{' num2str(p) '})']);
        axes('units','pixels','position',pos(axs1(p),:)),
        bar(means{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[means{p}-ses{p};means{p}+ses{p}],'color','k','linewidth',2)
        ylabel([{'SNR Model:'};{'ICMS Signal'};{'z-score'}],'fontsize',10,'fontname','arial')
        title(['Subject ' num2str(p)],'fontweight','bold','fontname','arial')
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
       
        axes('units','pixels','position',pos(axs2(p),:)),
        plot(beh_sum(1,:),'ob','markerfacecolor','b','markersize',5),hold on,
        plot(beh_sum(2,:),'sg','markerfacecolor','g','markersize',5)
        plot(beh_sum(3,:),'vr','markerfacecolor','r','markersize',5)
        plot(beh_sum(3,:),'r','linewidth',2)
        plot(beh_sum(2,:),'g','linewidth',2)
        plot(beh_sum(1,:),'b','linewidth',2)
        legend('Correct','False Alarm','Miss','Location',[0.715 0.124 .05 .05],'orientation','horizontal')
        
        ylim([0 100])
        xlim([0.5 size(beh_sum,2)+0.5])
        ylabel('Behavior (%)','fontsize',10,'fontname','arial')
        set(gca,'xtick',[1 size(beh_sum,2)])
        set(gca,'xticklabel',{'1';num2str(size(beh_sum,2))})
        set(gca,'ytick',[50 100])
        set(gca,'yticklabel',{'50';'100'})
        set(gca,'fontname','arial','fontsize',10)
        axis square
        
        axes('units','pixels','position',pos(axs3(p),:)),
        bar(meansC{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[meansC{p}-sesC{p};meansC{p}+sesC{p}],'color','k','linewidth',2)
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        ylabel([{'SNR Model:'};{'Internal Noise'};{'t-value'}],'fontsize',10,'fontname','arial')
        axis tight
        temp = get(gca,'ylim');
        ylim([0 temp(2)])
       
    elseif p == 3
        metric = eval(['length(file_list{' num2str(p) '})']);
        axes('units','pixels','position',pos(axs1(p),:)),
        bar(means{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[means{p}-ses{p};means{p}+ses{p}],'color','k','linewidth',2)
        title(['Subject ' num2str(p)],'fontweight','bold','fontname','arial')
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
       
        axes('units','pixels','position',pos(axs2(p),:)),
        plot(beh_sum(3,:),'vr','markerfacecolor','r','markersize',5),hold on,
        plot(beh_sum(2,:),'sg','markerfacecolor','g','markersize',5)
        plot(beh_sum(1,:),'ob','markerfacecolor','b','markersize',5)
        plot(beh_sum(3,:),'r','linewidth',2)
        plot(beh_sum(2,:),'g','linewidth',2)
        plot(beh_sum(1,:),'b','linewidth',2)

        ylim([0 100])
        xlim([0.5 size(beh_sum,2)+0.5])
        set(gca,'xtick',[1 size(beh_sum,2)])
        set(gca,'xticklabel',{'1';num2str(size(beh_sum,2))})
        set(gca,'ytick',[50 100])
        set(gca,'yticklabel',{'50';'100'})
        set(gca,'fontname','arial','fontsize',10)
        axis square
        
        axes('units','pixels','position',pos(axs3(p),:)),
        bar(meansC{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[meansC{p}-sesC{p};meansC{p}+sesC{p}],'color','k','linewidth',2)
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        xlabel('Session #','fontname','arial')
        axis tight
        temp = get(gca,'ylim');
        ylim([0 temp(2)])
        
    else
        metric = eval(['length(file_list{' num2str(p) '})']);
        axes('units','pixels','position',pos(axs1(p),:)),
        bar(means{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[means{p}-ses{p};means{p}+ses{p}],'color','k','linewidth',2)
        title(['Subject ' num2str(p)],'fontweight','bold','fontname','arial')
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
       
        axes('units','pixels','position',pos(axs2(p),:)),
        plot(beh_sum(3,:),'vr','markerfacecolor','r','markersize',5),hold on,
        plot(beh_sum(2,:),'sg','markerfacecolor','g','markersize',5)
        plot(beh_sum(1,:),'ob','markerfacecolor','b','markersize',5)
        plot(beh_sum(3,:),'r','linewidth',2)
        plot(beh_sum(2,:),'g','linewidth',2)
        plot(beh_sum(1,:),'b','linewidth',2)

        ylim([0 100])
        xlim([0.5 size(beh_sum,2)+0.5])
        set(gca,'xtick',[1 size(beh_sum,2)])
        set(gca,'xticklabel',{'1';num2str(size(beh_sum,2))})
        set(gca,'ytick',[50 100])
        set(gca,'yticklabel',{'50';'100'})
        set(gca,'fontname','arial','fontsize',10)
        axis square
        
        axes('units','pixels','position',pos(axs3(p),:)),
        bar(meansC{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[meansC{p}-sesC{p};meansC{p}+sesC{p}],'color','k','linewidth',2)
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
        temp = get(gca,'ylim');
        ylim([0 temp(2)])
       
    end
end

annotation('textbox',[.35 0.54 .05 .05],'string','Signal-to-Noise Models v. Behavior','fontname','arial','fontsize',10,'fontweight','bold','edgecolor','none');

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 11])
