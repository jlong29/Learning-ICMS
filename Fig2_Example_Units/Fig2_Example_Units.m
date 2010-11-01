%% Figure3: A neural correlate of ICMS Learning Behavior
clear

% Set dots per inch
dpi = 96;
figure('units','pixels','Position',[1 1 8.5 11].*dpi)

% From left to right
H = [1 2.5 4 5.5];

% From bottom to top
V = [4.5 5.75 7.0];

% Layout of axes
pos      = [H(1)  V(3)  1 1;
            H(2)  V(3)  1 1;
            H(3)  V(3)  1 1;
            H(4)  V(3)  1 1;
            
            H(1)  V(2)  1 1;
            H(2)  V(2)  1 1;
            H(3)  V(2)  1 1;
            H(4)  V(2)  1 1;
            
            H(1)  V(1)  1 1;
            H(2)  V(1)  1 1;
            H(3)  V(1)  1 1;
            H(4)  V(1)  1 1].*dpi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example Units for Stages 1, 2, and 3: axes 1-12 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic parameters
%%%%%%%%%%%%%%%%%%
samples = -200:32:800;  % us time axis for waveforms
edge    = 0:100;        % binning for ISI calculation (msec)
win     = 360;          % window for calculating firing rate in seconds
bin     = 10;           % bin for above window in seconds
r_isi   = 50;          % isi range for comparing distribution (msec)

% Load Unit data for all stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 1 Unit %
%%%%%%%%%%%%%%%%
unit1 = 'S1_V1_sig003a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V1_092408_1HD_ICMS_tr_clean_Ex.mat',unit1)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
load('V1_092408_spikes&Events_1_16.mat',unit1(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes1 = ' sp_list{1} ';'])
eval(['wf1     = ' wvf_list{1} ';'])

% scale from V to mV
wf1 = wf1.*1000;
clear wvf_list* sp_list* S* sig*

% Stage 2 Unit %
%%%%%%%%%%%%%%%%
unit2 = 'S2_V7_sig007a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V7_032909_baseline_and_ICMS_clean_Ex.mat',unit2)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
load('V7_032909_spikes&Events1_16.mat',unit2(end-6:end))
wvf_list  = who('S*');
sp_list   = who('sig*');

eval(['spikes2 = ' sp_list{1} ';'])
eval(['wf2     = ' wvf_list{1} ';'])

% scale from V to mV
wf2 = wf2.*1000;

clear wvf_list* sp_list* S* sig*

% Stage 3 Unit %
%%%%%%%%%%%%%%%%
unit3 = 'S3_V1_sig002a';
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
load('V1_093008_ICMS_clean_Ex.mat',unit3)
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
load('V1_093008_spikes&Events_1_16.mat',unit3(end-6:end))
wvf_list = who('S*');
sp_list  = who('sig*');

eval(['spikes3 = ' sp_list{1} ';'])
eval(['wf3     = ' wvf_list{1} ';'])

% scale from V to mV
wf3 = wf3.*1000;

clear wvf_list* sp_list* S* sig*

%%%%%%%%%%%%%%%%%%%%
%%% Stage 1 Unit %%%
%%%%%%%%%%%%%%%%%%%%
start   = wf1(1:1000,:);
isi1    = 1000.*diff(spikes1(1:5001));

fin     = wf1(end-999:end,:);
isi2    = 1000.*diff(spikes1(end-5000:end));

meanst     = mean(start);
strss      = sum((start - repmat(meanst,1000,1)).^2,2);

meanfin    = mean(fin);
finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);

[c,stind]  = sort(strss);
[d,finind] = sort(finrss);

% Axes 1
% Look at early session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = axes('units','pixels','position',pos(1,:));
    plot(start(stind(1:200),:)','k')
    xlim([1 32])
    temp  = logical(unit1 == '_');
    unit1(temp) = ' ';
    ylabel([{'Stage 1 Unit'};{'uV'}])
    title('Early (E)','fontweight','bold')
    set(h,'xticklabel','')
    set(h,'ylim',[-200 200])
    y = get(h,'ylim');
    ticks  = get(h,'ytick');
    set(h,'ytick',[ticks(1) 0 ticks(end)])
    ticks  = get(h,'ytick');
    set(gca,'fontname','arial','fontsize',10)
    axis square

% Axes 2
% Look at late session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = axes('units','pixels','position',pos(2,:));
    plot(fin(finind(1:200),:)','color',[0.3912 0.3990 0.350])
    xlim([1 32])
    title('Late (L)','fontweight','bold')
    set(j,'xticklabel','')
    set(j,'ylim',y)
    set(j,'ytick',ticks)
    set(j,'yticklabel','')
    set(gca,'fontname','arial','fontsize',10)
    axis square

% Axes 3
% Look at early and late session ISIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use KS statistic for comparing isi distributions within range 1-100ms
bnd1  = logical(isi1 <= r_isi);
bnd2  = logical(isi2 <= r_isi);
axes('units','pixels','position',pos(3,:)),
               plot(edge,histc(isi1,edge),'k','linewidth',2),hold on,
               plot(edge,histc(isi2,edge),'color',[0.3912 0.3990 0.350],'linewidth',2)
               axis([-5 50 0 75])
               set(gca,'xtick',[0 50])
               set(gca,'xticklabel','0|50')
               ylabel('Counts')
               title('ISI Distribution','fontweight','bold')
               set(gca,'fontname','arial','fontsize',10)
               axis square
% Axes 4
% Look at early and late session firing rate of unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line1 = spikes1(1):bin:spikes1(1) + win;

unit_fr1   = zeros(1,length(time_line1)-1);
data1      = spikes1;

for p = 2:length(time_line1)

    ind = find(data1 >= time_line1(p-1) & data1 < time_line1(p));
    counts = length(ind);

    % Speed up
    unit_fr1(p-1) = unit_fr1(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr1   = unit_fr1 ./ bin;

% Look at early and late session firing rate of unit
time_line2 = spikes1(end) - win:bin:spikes1(end);

unit_fr2   = zeros(1,length(time_line2)-1);
data1      = spikes1;

for p = 2:length(time_line2)

    ind = find(data1 >= time_line2(p-1) & data1 < time_line2(p));
    counts = length(ind);

    % Speed up
    unit_fr2(p-1) = unit_fr2(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr2   = unit_fr2 ./ bin;

mean1      = mean(unit_fr1);
std1       = std(unit_fr1);

mean2      = mean(unit_fr2);
std2       = std(unit_fr2);

axes('units','pixels','position',pos(4,:))
    bar(1,mean1,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),hold on,
    bar(2,mean2,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),
    line([1 2;1 2],[mean1-std1 mean2-std2;mean1+std1 mean2+std2],'color','k','linewidth',3)
    xlim([0.5 2.5])
    ylabel('spikes/s')
    title('Firing Rate','fontweight','bold')
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'E';'L'})
    set(gca,'fontname','arial','fontsize',10)
    axis square

%%%%%%%%%%%%%%%%%%%%
%%% Stage 2 Unit %%%
%%%%%%%%%%%%%%%%%%%%
start   = wf2(1:1000,:);
isi1    = 1000.*diff(spikes2(1:5001));

fin     = wf2(end-999:end,:);
isi2    = 1000.*diff(spikes2(end-5000:end));

meanst     = mean(start);
strss      = sum((start - repmat(meanst,1000,1)).^2,2);

meanfin    = mean(fin);
finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);

[c,stind]  = sort(strss);
[d,finind] = sort(finrss);

% Axes 5
% Look at early session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = axes('units','pixels','position',pos(5,:));
    plot(start(stind(1:200),:)','k')
    xlim([1 32])
    temp   = logical(unit2 == '_');
    label  = unit2;
    label(temp) = ' ';
    ylabel([{'Stage 2 Unit'};{'uV'}])
    set(h,'xticklabel','')
    y      = get(h,'ylim');
    set(h,'ylim',y)
    ticks  = get(h,'ytick');
    set(h,'ytick',[ticks(1) 0 ticks(end)])
    ticks  = get(h,'ytick');
    set(gca,'fontname','arial','fontsize',10)
    axis square

% Axes 6
% Look at late session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = axes('units','pixels','position',pos(6,:));
    plot(fin(finind(1:200),:)','color',[0.3912 0.3990 0.350])
    xlim([1 32])
    set(j,'xticklabel','')
    set(j,'ylim',y)
    set(j,'ytick',ticks)
    set(j,'yticklabel','')
    axis square

% Axes 7
% Look at early and late session ISIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd1  = logical(isi1 <= r_isi);
bnd2  = logical(isi2 <= r_isi);
[h,p] = kstest2(isi1(bnd1),isi2(bnd2),.01,'unequal');
axes('units','pixels','position',pos(7,:)),
               plot(edge,histc(isi1,edge),'k','linewidth',2),hold on,
               plot(edge,histc(isi2,edge),'color',[0.3912 0.3990 0.350],'linewidth',2)
               axis([-5 50 0 75])
               ylabel('Counts')
               set(gca,'xtick',[0 50])
               set(gca,'xticklabel','0|50')
               set(gca,'fontname','arial','fontsize',10)
               axis square
% Axes 8
% Look at early and late session firing rate of unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line1 = spikes2(1):bin:spikes2(1) + win;

unit_fr1   = zeros(1,length(time_line1)-1);
data1      = spikes2;

for p = 2:length(time_line1)

    ind = find(data1 >= time_line1(p-1) & data1 < time_line1(p));
    counts = length(ind);

    % Speed up
    unit_fr1(p-1) = unit_fr1(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr1   = unit_fr1 ./ bin;

% Look at early and late session firing rate of unit
if strcmpi(label,'S2 R29 sig022a')
    time_line2 = 2200:bin:2560;
else
    time_line2 = spikes2(end) - win:bin:spikes2(end);
end

unit_fr2   = zeros(1,length(time_line2)-1);
data1      = spikes2;

for p = 2:length(time_line2)

    ind = find(data1 >= time_line2(p-1) & data1 < time_line2(p));
    counts = length(ind);

    % Speed up
    unit_fr2(p-1) = unit_fr2(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr2   = unit_fr2 ./ bin;

mean1      = mean(unit_fr1);
std1       = std(unit_fr1);

mean2      = mean(unit_fr2);
std2       = std(unit_fr2);

axes('units','pixels','position',pos(8,:))
    bar(1,mean1,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),hold on,
    bar(2,mean2,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),
    line([1 2;1 2],[mean1-std1 mean2-std2;mean1+std1 mean2+std2],'color','k','linewidth',3)
    xlim([0.5 2.5])
    ylabel('spikes/s')
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'E';'L'})
    set(gca,'fontname','arial','fontsize',10)
    axis square


%%%%%%%%%%%%%%%%%%%%
%%% Stage 3 Unit %%%
%%%%%%%%%%%%%%%%%%%%
start   = wf3(1:1000,:);
isi1    = 1000.*diff(spikes3(1:5001));

fin     = wf3(end-999:end,:);
isi2    = 1000.*diff(spikes3(end-5000:end));

meanst     = mean(start);
strss      = sum((start - repmat(meanst,1000,1)).^2,2);

meanfin    = mean(fin);
finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);

[c,stind]  = sort(strss);
[d,finind] = sort(finrss);

% Axes 9
% Look at early session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = axes('units','pixels','position',pos(9,:));
    plot(start(stind(1:200),:)','k')
    xlim([1 32])
    temp   = logical(unit3 == '_');
    unit3(temp) = ' ';
    xlabel('Time (\musec)')
    ylabel([{'Stage 3 Unit'};{'uV'}])
    set(h,'xtick',[1 32])
    set(h,'xticklabel','-200|600')
    set(h,'ylim',[-200 200])
    set(h,'ytick',[-200 0 200])
    set(h,'yticklabel','-200|0|200')
    ticks  = get(h,'ytick');
    lticks = get(h,'yticklabel');
    set(gca,'fontname','arial','fontsize',10)
    axis square

% Axes 10
% Look at late session waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = axes('units','pixels','position',pos(10,:));
    plot(fin(finind(1:200),:)','color',[0.3912 0.3990 0.350])
    xlim([1 32])
    xlabel('Time (\musec)')
    set(j,'xtick',[1 32])
    set(j,'xticklabel','-200|600')
    set(j,'ylim',y)
    set(j,'ytick',ticks)
    set(j,'yticklabel','')
    axis square

% Axes 11
% Look at early and late session ISIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd1  = logical(isi1 <= r_isi);
bnd2  = logical(isi2 <= r_isi);
[h,p] = kstest2(isi1(bnd1),isi2(bnd2),.01,'unequal');
axes('units','pixels','position',pos(11,:)),
               plot(edge,histc(isi1,edge),'k','linewidth',2),hold on,
               plot(edge,histc(isi2,edge),'color',[0.3912 0.3990 0.350],'linewidth',2)
               axis([-5 50 0 120])
               xlabel('Time (msec)')
               ylabel('Counts')
               set(gca,'xtick',[0 50])
               set(gca,'xticklabel','0|50')
               set(gca,'fontname','arial','fontsize',10)
               axis square
% Axes 12
% Look at early and late session firing rate of unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_line1 = spikes3(1):bin:spikes3(1) + win;

unit_fr1   = zeros(1,length(time_line1)-1);
data1      = spikes3;

for p = 2:length(time_line1)

    ind = find(data1 >= time_line1(p-1) & data1 < time_line1(p));
    counts = length(ind);

    % Speed up
    unit_fr1(p-1) = unit_fr1(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr1   = unit_fr1 ./ bin;

% Look at early and late session firing rate of unit
time_line2 = spikes3(end) - win:bin:spikes3(end);

unit_fr2   = zeros(1,length(time_line2)-1);
data1      = spikes3;

for p = 2:length(time_line2)

    ind = find(data1 >= time_line2(p-1) & data1 < time_line2(p));
    counts = length(ind);

    % Speed up
    unit_fr2(p-1) = unit_fr2(p-1) + counts;
    if ~isempty(ind)
        data1    = data1(ind(end):end);
    end
end

% Convert to spike per second
unit_fr2   = unit_fr2 ./ bin;

mean1      = mean(unit_fr1);
std1       = std(unit_fr1);

mean2      = mean(unit_fr2);
std2       = std(unit_fr2);

axes('units','pixels','position',pos(12,:))
    bar(1,mean1,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),hold on,
    bar(2,mean2,'facecolor',[0.3912 0.3990 0.350],'edgecolor','k'),
    line([1 2;1 2],[mean1-std1 mean2-std2;mean1+std1 mean2+std2],'color','k','linewidth',3)
    xlim([0.5 2.5])
    ylabel('spikes/s')
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'E';'L'})
    set(gca,'fontname','arial','fontsize',10)
    axis square

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')