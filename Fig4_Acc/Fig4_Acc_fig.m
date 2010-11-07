% Fig4: Accelerometer Data: obvert behavior
close all
clear

%%%%%%%%%%%%%%%%%%%%%
%%% Figure layout %%%
%%%%%%%%%%%%%%%%%%%%%
% Set dots per inch
dpi = 96;
figure('units','pixels','Position',[1 1 8.5 7.6].*dpi)

% From left to right
H = 2.5:1.5:5+1.5;

% From top to bottom
V = 2.25*2+0.5:-1.50:1.0;

% Layout of axes
pos      = [H(1)  V(1)  1.0 1.0
            H(1)  V(2)  1.0 1.0
            H(1)  V(3)  1.0 1.0
            
            H(2)  V(1)  1.0 1.0
            H(2)  V(2)  1.0 1.0
            H(2)  V(3)  1.0 1.0
            
            H(3)  V(1)  1.0 1.0
            H(3)  V(2)  1.0 1.0
            H(3)  V(3)  1.0 1.0].*dpi;

% Load data
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig4_Acc')
load('Fig4_Acc_figure_data')

% Subject IDs
subjects = {'R29*','V4*','V8*'};
win   = 1000;
time  = -win:win;

% Grab the data
vars1 = who(subjects{1});
vars2 = who(subjects{2});
vars3 = who(subjects{3});

for q = 1:length(subjects)
    % Grab the data
    data = who(subjects{q});

    for n = 1:length(data)/2
        
        accm   = eval(data{2*n - 1});
        accm   = accm(3,:);
        
        accse  = eval(data{2*n});
        accse  = accse(3,:);
        
        axes('units','pixels','position',pos(3*(q-1) +n,:),'fontname','arial','fontsize',10)
        plot(time,accm,'k','linewidth',2),hold on
        ciplot(accm-accse,accm+accse,time,'k',0.3)
        set(gca,'xtick',[-1000 -500 0 500 1000],'xticklabel','-1000||0||1000')
        set(gca,'ticklength',[.025 .1])
        
        % get y limits and set title
        if n == 1
            ymax = mean(accm) + 5*mean(accse);
            ymin = mean(accm) - 2.5*mean(accse);
            
            title([{['Subject' num2str(q) ':']};{'response to ICMS'}],'fontname','arial','fontsize',10)
        end
        
        % set y limits
        set(gca,'ylim',[ymin ymax])
        
        % set xlabel
        if n == 3
            xlabel('Time (msec)','fontname','arial','fontsize',10)
        end
        
        if q == 1
            ylabel([{'Acceleration'};{'z-score'}],'fontname','arial','fontsize',10)
        end
    end
end

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 11])
