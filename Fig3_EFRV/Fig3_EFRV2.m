%% Figure3 cont: A neural correlate of ICMS Learning Behavior
clear

scrnsz = get(0,'screensize');
figure('Position',scrnsz,'units','pixels')

% Set dots per inch
dpi = 96;

% From left to right
H = [3.5 5 6.5 8 9.5];

% From bottom to top
V = [4 5.25];

% Layout of axes
pos      = [H(1)*dpi  V(2)*dpi  1*dpi 1*dpi;
            H(2)*dpi  V(2)*dpi  1*dpi 1*dpi;
            H(3)*dpi  V(2)*dpi  1*dpi 1*dpi;
            H(4)*dpi  V(2)*dpi  1*dpi 1*dpi;
            H(5)*dpi  V(2)*dpi  1*dpi 1*dpi;
            
            H(1)*dpi  V(1)*dpi  1*dpi 1*dpi;
            H(2)*dpi  V(1)*dpi  1*dpi 1*dpi;
            H(3)*dpi  V(1)*dpi  1*dpi 1*dpi;
            H(4)*dpi  V(1)*dpi  1*dpi 1*dpi;
            H(5)*dpi  V(1)*dpi  1*dpi 1*dpi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ensemble NC Learning Subplots: Axes 1 - 10 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
% Grab ensemble Firing rate data for test hemisphere
load('Subject_SDs_test.mat')
% Grab subjects' behavioral data
load('Subjects_Behavior.mat')

means = {R29_mean; V1_mean; V4_mean; V7_mean; V8_mean};
ses   = {R29_SE;   V1_SE;   V4_SE;   V7_SE;   V8_SE};

% Axes pointers
axs1   = 1:5;
axs2   = 6:10;

for p = 1:size(pos,1)/2
    beh_sum = subject_Behavior{p};
    
    if p == 1
        metric = eval(['length(file_list' num2str(p) ')']);
        axes('units','pixels','position',pos(axs1(p),:)),
        bar(means{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[means{p}-ses{p};means{p}+ses{p}],'color','k','linewidth',2)
        ylabel([{'Ensemble Firing'};{'Rate Variablity'}],'fontsize',10,'fontname','arial')
        title(['Subject ' num2str(p)],'fontweight','bold','fontname','arial')
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
%         ylim([0 1])
%         set(gca,'ytick',[0 1])
        
        axes('units','pixels','position',pos(axs2(p),:)),
        plot(beh_sum(3,:),'vr','markerfacecolor','r','markersize',3),hold on,
        plot(beh_sum(2,:),'sg','markerfacecolor','g','markersize',3)
        plot(beh_sum(1,:),'ob','markerfacecolor','b','markersize',3)
        plot(beh_sum(3,:),'r','linewidth',2)
        plot(beh_sum(2,:),'g','linewidth',2)
        plot(beh_sum(1,:),'b','linewidth',2)

        ylim([0 100])
        xlim([0.5 size(beh_sum,2)+0.5])
        ylabel([{'Relative'}; {'Behavior (%)'}],'fontsize',10,'fontname','arial')
        set(gca,'xtick',[1 size(beh_sum,2)])
        set(gca,'xticklabel',{'1';num2str(size(beh_sum,2))})
        set(gca,'ytick',[50 100])
        set(gca,'yticklabel',{'50';'100'})
        set(gca,'fontname','arial','fontsize',10)
        axis square
    else
        metric = eval(['length(file_list' num2str(p) ')']);
        axes('units','pixels','position',pos(axs1(p),:)),
        bar(means{p},'edge','k','facecolor',[0.3912 0.3990 0.350],'barwidth',1)
        line([1:metric;1:metric],[means{p}-ses{p};means{p}+ses{p}],'color','k','linewidth',2)
        title(['Subject ' num2str(p)],'fontweight','bold','fontname','arial')
        set(gca,'xtick',[1 length(means{p})])
        set(gca,'xticklabel',{'1';num2str(metric)})
        axis tight
%         ylim([0 1])
%         set(gca,'ytick',[0 1])
        
        axes('units','pixels','position',pos(axs2(p),:)),
        plot(beh_sum(3,:),'vr','markerfacecolor','r','markersize',3),hold on,
        plot(beh_sum(2,:),'sg','markerfacecolor','g','markersize',3)
        plot(beh_sum(1,:),'ob','markerfacecolor','b','markersize',3)
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
    end
end

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')

annotation('textbox',[.405 0.85 .05 .05],'string','Ensemble Firing Rate Variability v. Behavior','fontname','arial','fontsize',10,'fontweight','bold','edgecolor','none');
annotation('textbox',[.497 0.45 .05 .05],'string','Session #', 'fontname','arial','fontsize',10,'edgecolor','none');