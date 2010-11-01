%%% Figure4_final: Ensemble Firing Rate Variability v. Accelerometer %%%

dpi   = 96;
figure('units','pixels','Position',[1 1 8.5*dpi 11*dpi])

H     = [1 4.25];
V     = 4;
pos   = [H(1)*dpi   V(1)*dpi 2.5*dpi 2.58*dpi;
         H(2)*dpi   V(1)*dpi 3*dpi 2.58*dpi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Accelerometer and Spike Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('c:\program files\matlab\r2006b\work\Learning_ICMS_figures')
load('Acc_and_Spike_example.mat')

% Lowpass filter Acc data at 50 Hz
[B,A] = butter(5,50/1000,'low');
Acc= filter(B,A,Acc);

axes('units','pixels','position',pos(1,:)),
plot(time_line,HF_Spikes(:,1)+10,'k'),hold on
plot(time_line,HF_Spikes(:,4)+8,'k')
plot(time_line,HF_Spikes(:,5)+6,'k')

plot(time_line2,Acc(:,1)+4.25,'b','linewidth',2)
plot(time_line2,Acc(:,3)+2,'b','linewidth',2)
plot(time_line2,Acc(:,2)-.6,'b','linewidth',2)
title([{'Simultaneous Spike and'}; {'Accelerometer Data'}],'fontsize',12,'fontname','arial','fontweight','bold')
ylabel('Arbitrary Units','fontsize',12,'fontname','arial')
xlabel('Time (msec)','fontsize',12,'fontname','arial')
set(gca,'xtick',[0 length(time_line)/2 length(time_line)],'xticklabel','0|500|1000')
axis tight
temp = get(gca,'ylim');
set(gca,'xlim',[0 length(time_line)],'ylim',[temp(1)-0.1 temp(2)])
set(gca,'yticklabel','')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Accelerometer Null Hypothesis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
load('Null_Acc2EfrV_2.mat')
thres = 5;

% R^2 of Ensemble Firing Rate Variability v. Acceleration events
x = [subject_Acc{1,1}(thres,:)',subject_Acc{1,2}';
     subject_Acc{2,1}(thres,:)',subject_Acc{2,2}';
     subject_Acc{3,1}(thres,:)',subject_Acc{3,2}';
     subject_Acc{4,1}(thres,:)',subject_Acc{4,2}';
     subject_Acc{5,1}(thres,:)',subject_Acc{5,2}'];

R2 = corr(x).^2;
R2 = R2(2);

axes('units','pixels','position',pos(2,:)),
% Each Subject with different Marker and Mark Active Exploration Sessions
plot(subject_Acc{1,1}(thres,:),subject_Acc{1,2},'ok','markerfacecolor','k','markersize',5),hold on
plot(subject_Acc{2,1}(thres,:),subject_Acc{2,2},'sk','markerfacecolor','k','markersize',5)
plot(subject_Acc{3,1}(thres,:),subject_Acc{3,2},'dk','markerfacecolor','k','markersize',5)
plot(subject_Acc{4,1}(thres,:),subject_Acc{4,2},'^k','markerfacecolor','k','markersize',5)
plot(subject_Acc{5,1}(thres,:),subject_Acc{5,2},'<k','markerfacecolor','k','markersize',5)

plot(subject_Acc{1,1}(thres,6),subject_Acc{1,2}(6),'ob','markerfacecolor','c','markersize',5)
plot(subject_Acc{2,1}(thres,3),subject_Acc{2,2}(3),'sb','markerfacecolor','c','markersize',5)
plot(subject_Acc{3,1}(thres,4),subject_Acc{3,2}(4),'db','markerfacecolor','c','markersize',5)
plot(subject_Acc{4,1}(thres,5),subject_Acc{4,2}(5),'^b','markerfacecolor','c','markersize',5)
plot(subject_Acc{5,1}(thres,4),subject_Acc{5,2}(4),'<b','markerfacecolor','c','markersize',5)
title([{'Ensemble Firing Rate Variablity'};{'vs. Acceleration Events'}],'fontsize',12,'fontweight','bold','fontname','arial')
ylabel('Ensemble Firing Rate Variability','fontname','arial','fontsize',10)
xlabel('Acceleration Events','fontname','arial','fontsize',10)
axis tight

temp1 = get(gca,'ylim');
set(gca,'ylim',[temp1(1)-0.1 temp1(2)+0.1])
temp2 = get(gca,'xlim');
set(gca,'xlim',[temp2(1)-100 temp2(2)+200])
set(gca,'ytick',[3 5 7])
set(gca,'xtick',[0 900 1800],'xticklabel','0|900|1800')

text(1300, 3.5,['R^2 = ' num2str(R2,1)],'fontname','arial','fontsize',10,'fontweight','bold')

set(findobj(gcf,'Type','axes'),'Fontname','arial')
set(findobj(gcf,'Type','axes'),'box','off')
h = legend('S1','S2','S3','S4','S5','location','northeast');
set(h,'edgecolor','k')
set(gcf,'color','w')
