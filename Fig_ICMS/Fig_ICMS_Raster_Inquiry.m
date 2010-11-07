%%%%%%%
% R29 %
%%%%%%%
disp('V8_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list1] = find_files('V8','1_16');

% Generate SNR signal and internal noise models
for n = 1:length(file_list1)
    % Let's grab the file date
    file  = file_list1{n};
    [temp,name] = fileparts(file);
    parse = findstr(name,'_');
    date  = name(parse(1)+1:parse(2)-1);
    
    % Generate PSTHs for each unit reference to ICMS
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file,'Event015',10,500,500);
    
    [all_rasters,ave_fr,unitIds,event_time] = units_all_rasters(file,'Event015',500,500);
    
%     [all_rasters,ave_fr,unitIds,event_time] = units_all_rasters(file,'Event015',100,100);
    pause
    close all
    
end

% set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 11])
% print -dpng -r600 test

%%
subject = 'V8';
cd(['C:\Data\Plexon_Data\Test_Data\Learning\' subject])
[file_list1] = find_files(subject,'1_16');

scrnsz = get(0,'screensize');
figure('Position',scrnsz)
for n = 1:length(file_list1)
    % Let's grab the file date
    file  = file_list1{n};
    [temp,name] = fileparts(file);
    parse = findstr(name,'_');
    date  = name(parse(1)+1:parse(2)-1);
    
    % Look at paired differences between start trial samples and start
    % trial/pre-stim sample
    [pd_base,pd_stim,unitIds] = change_pre_stim_fr(file,1000);
    
    upper = max([pd_stim(:);pd_base(:)]);
    lower = min([pd_stim(:);pd_base(:)]);
    edge = lower:upper;
    
    subplot(211),
    hist(pd_base(:),edge),
    hold on,hist(pd_stim(:),edge),
    kids = get(gca,'children');
    set(kids(2),'Edgecolor','k','Facecolor','k','EdgeAlpha',1,'FaceAlpha',1)
    set(kids(1),'Edgecolor',[0.3912 0.3990 0.350],'Facecolor',...
        [0.3912 0.3990 0.350],'EdgeAlpha',0.8,'FaceAlpha',0.8)
    l = legend('\Deltabase','\Deltastim','location','northwest');
    xlim([-40 40])
    
    title([{[subject ' ' date]};{['Rank Sum p-value = ' num2str(ranksum(pd_base(:),pd_stim(:)))]}])
    ylabel('Counts','fontname','arial','fontsize',10)
    xlabel('Data Limits','fontname','arial','fontsize',10)
    axis square
    
    temp1 = pd_base(:);
    temp2 = pd_stim(:);
    [temp,ind]=sort(temp1);
    subplot(212),
    plot(temp1(ind),temp2(ind),'.k','markersize',5)
    ylim([-40 40])
    xlim([-40 40])
    
    p1 = polyfit(temp1,temp2,0);
    p2 = polyfit(temp2,temp1,0);
    line(get(gca,'xlim'),[p1 p1],'color','k','linestyle','--')
    line([p2 p2],get(gca,'ylim'),'color','k','linestyle','--')
    
    title('\Deltabase vs \Deltastim','fontname','arial','fontsize',10)
    xlabel('\Deltabase data limits','fontname','arial','fontsize',10)
    ylabel('\Deltastim data limits','fontname','arial','fontsize',10)
    axis square
    pause
    clf
end
