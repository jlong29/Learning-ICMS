%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure ICMS examples %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear

%%%%%%%%%%%%%%%%%%%%%
%%% Figure layout %%%
%%%%%%%%%%%%%%%%%%%%%
% Set dots per inch
dpi = 96;
figure('units','pixels','Position',[1 1 8.5 7.6].*dpi)

% From left to right
H = [1 2.10 3.20 4.30 5.40 6.5 3.30 5.50];

% From bottom to top
V = [6.3 5.0 2.75 0.75 2.20];

% Layout of axes
pos      = [H(1)  V(1)  0.75 0.75
            H(1)  V(2)  0.75 0.75
            
            H(2)  V(1)  0.75 0.75
            H(2)  V(2)  0.75 0.75
            
            H(3)  V(1)  0.75 0.75
            H(3)  V(2)  0.75 0.75
            
            H(4)  V(1)  0.75 0.75
            H(4)  V(2)  0.75 0.75
            
            H(5)  V(1)  0.75 0.75
            H(5)  V(2)  0.75 0.75
            
            H(6)  V(1)  0.75 0.75
            H(6)  V(2)  0.75 0.75
            
            H(1)  V(3)  1.5  1.5
            H(1)  V(4)  1.5  1.5
            
            H(7)  V(5)  1.75  1.75
            H(8)  V(5)  1.75  1.75].*dpi;

%%%%%%%%%%%%%%%%%%%%%
%%% Example Units %%%
%%%%%%%%%%%%%%%%%%%%%
%%% IMPORT DATA %%%
% Unit1
unit1 = 'S2_R29_sig017b';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('R29_051008_rec1HD_blockICMS_clean_Ex.mat',unit1)
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
load('R29_051008_spike&Events17_32.mat',unit1(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes1 = ' sp_list{1} ';'])
eval(['wf1     = ' wvf_list{1} ';'])

% scale from V to mV
wf1 = wf1.*1000;
clear wvf_list* sp_list* S* sig*

% Unit2
unit2 = 'S2_V7_sig008a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V7_032909_baseline_and_ICMS_clean_Ex.mat',unit2)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
load('V7_032909_spikes&Events1_16.mat',unit2(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes2 = ' sp_list{1} ';'])
eval(['wf2     = ' wvf_list{1} ';'])

% scale from V to mV
wf2 = wf2.*1000;
clear wvf_list* sp_list* S* sig*

% Unit3
unit3 = 'S2_V7_sig016a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V7_032909_baseline_and_ICMS_clean_Ex.mat',unit3)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
load('V7_032909_spikes&Events1_16.mat',unit3(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes3 = ' sp_list{1} ';'])
eval(['wf3     = ' wvf_list{1} ';'])

% scale from V to mV
wf3 = wf3.*1000;
clear wvf_list* sp_list* S* sig*

% Unit4
unit4 = 'S2_V8_sig005a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V8_033109_baseline_and_ICMS_clean_Ex.mat',unit4)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
load('V8_033109_spikes&Events1_16.mat',unit4(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes4 = ' sp_list{1} ';'])
eval(['wf4     = ' wvf_list{1} ';'])

% scale from V to mV
wf4 = wf4.*1000;
clear wvf_list* sp_list* S* sig*

% Unit5
unit5 = 'S2_V8_sig011a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V8_033109_baseline_and_ICMS_clean_Ex.mat',unit5)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
load('V8_033109_spikes&Events1_16.mat',unit5(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes5 = ' sp_list{1} ';'])
eval(['wf5     = ' wvf_list{1} ';'])

% scale from V to mV
wf5 = wf5.*1000;
clear wvf_list* sp_list* S* sig*

% Unit6
unit6 = 'S2_V8_sig013a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V8_033109_baseline_and_ICMS_clean_Ex.mat',unit6)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
load('V8_033109_spikes&Events1_16.mat',unit6(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes6 = ' sp_list{1} ';'])
eval(['wf6     = ' wvf_list{1} ';'])

% scale from V to mV
wf6 = wf6.*1000;
clear wvf_list* sp_list* S* sig*

%%% Consolidate data
all_spikes = {spikes1,spikes2,spikes3,spikes4,spikes5,spikes6};
all_wfs    = {wf1,wf2,wf3,wf4,wf5,wf6};

clear spikes* wf* unit*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESS FIGURE %%%

for n = 1:length(all_spikes)
    % Select data
    data   = all_wfs{n}(1:1000,:);
    isi    = 1000.*diff(all_spikes{n}(1:5001));

    %%%%%%%%%%%%%%%%%
    %%% Axes 1    %%%
    %%% waveforms %%%
    h = axes('units','pixels','position',pos(2*(n-1)+1,:));
        plot(data(1:200,:)','k')
        xlim([1 32])
        set(h,'xticklabel','')
        set(h,'ylim',[-250 250])
        y = get(h,'ylim');
        set(h,'ytick',[-200 0 200])
        set(gca,'fontname','arial','fontsize',8)
        set(h,'xtick',[1 32])
        set(h,'xticklabel','-200|600')
        title('Waveforms','fontsize',10,'fontname','arial','fontweight','bold')
        xlabel('Time (\musec)','fontname','arial','fontsize',8)
        
        if n == 1
            ylabel('\muV','fontname','arial','fontsize',10)
        end
        
    %%%%%%%%%%%%%%
    %%% Axes 2 %%%
    %%% ISIs   %%%
    edge    = 0:100;        % binning for ISI calculation (msec)
    r_isi   = 50;          % isi range for comparing distribution (msec)
    bnd  = logical(isi <= r_isi);

    axes('units','pixels','position',pos(2*(n-1)+2,:)),
       h = histc(isi,edge);
       plot(edge,h,'k','linewidth',1),
       xlim([-5 30])
       ylim([0 max(h)+5-mod(max(h),5)+5])
       line([0 0],[0 max(h)+5-mod(max(h),5)+5],'color',[0.3912 0.3990 0.350],...
           'linestyle','--')
       set(gca,'xtick',[0 25],'ytick',[0 max(h)+5-mod(max(h),5)])
       set(gca,'xticklabel','0|25')
       xlabel('ISI (ms)','fontname','arial','fontsize',8)
       set(gca,'fontname','arial','fontsize',8)
       
       if n == 1
           ylabel('Counts','fontname','arial','fontsize',10)
       end
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ICMS Response Example %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load example file
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig_ICMS')
load('Fig_ICMS_example')

% raster template from 'units_all_rasters.m'
ind     = 1:1001;
firings =[];
for t = ind
    fired=find(raster_data(t,:))';
    firings=[firings; t+0*fired,fired];
end

axes('units','pixels','position',pos(13,:)),
plot(firings(:,1),firings(:,2),'.k','markersize',1),hold on
line([501 501],[0 size(raster_data,2)+1],'color','r','linestyle','--')

title('ICMS Response: Raster','fontname','arial','fontsize',10),
ylabel('Trials','fontname','arial','fontsize',10)

xlim([0 size(raster_data,1)]),ylim([0 size(raster_data,2)])
set(gca,'yticklabel','')
set(gca,'xtick',[1 251 501 751 1001])
set(gca,'xticklabel','-500|-250|0|250|500')
set(gca,'fontname','arial','fontsize',10)
axis square
    
% PSTH template from 'units_all_PSTH.m'
ave_fr     = mean(psth_data(1:50));
event_time = -500:10:490;

axes('units','pixels','position',pos(14,:)),
bar(event_time,psth_data,'facecolor','k','edgecolor','k'),hold on
xlim([-500 490])
set(gca,'xtick',[-500 -250 -5 250 490])
set(gca,'xticklabel','-500|-250|0|250|500')
set(gca,'fontsize',8)

title('ICMS Response: PSTH','fontname','arial','fontsize',10),
ylabel([{'Rate'};{'(Spikes/second)'}],'fontname','arial','fontsize',10)
xlabel('Time (ms)','fontname','arial','fontsize',10)
axis square    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR Models: Increasing Signal and Decreasing Noise Axes 1 - 4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the relevant files
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig_ICMS')
file = 'Figure_ICMS_PSTH_SNR_bin10_win500.mat';
load(file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of Unit excited by ICMS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line = -500:10:490;
axes('units','pixels','position',pos(15,:)),
    bar(time_line,R29_ICMS_psth_2(:,12),'k','barwidth',1)
    xlim([-500 500])
    title([{'Signal:'}; {'Excitation'}],'fontname','arial','fontsize',10,'fontweight','bold')
    ylabel([{'Rate'};{'(Spikes/second)'}],'fontname','arial','fontsize',10)
    xlabel('Time (ms)','fontname','arial','fontsize',10)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of Unit inhibited by ICMS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line = -500:10:490;
axes('units','pixels','position',pos(16,:)),
    bar(time_line,V7_ICMS_psth_1(:,11),'k','barwidth',1)
    xlim([-500 500])
    title([{'Internal Noise:'}; {'Inhibition'}],'fontname','arial','fontsize',10,'fontweight','bold')
    xlabel('Time (ms)','fontname','arial','fontsize',10)
    
annotation('textbox',[0.4900 0.560 .1 .05],'VerticalAlignment','top','string',...
    'Signal-to-Noise Ratio Models','linestyle','none',...
    'fontname','arial','fontsize',10,'fontweight','bold');


set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 7.6])

% print -dpng -r600 -opengl Figure_ICMS