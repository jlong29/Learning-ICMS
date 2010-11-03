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
subject = 'R29';
cd(['C:\Data\Plexon_Data\Test_Data\Learning\' subject])
[file_list1] = find_files(subject,'1_16');

% Generate SNR signal and internal noise models
ratio = cell(length(file_list1),1);
for n = 1:length(file_list1)
    % Let's grab the file date
    file  = file_list1{n};
    [temp,name] = fileparts(file);
    parse = findstr(name,'_');
    date  = name(parse(1)+1:parse(2)-1);
    
    % Generate PSTHs for each unit reference to ICMS
    [pre_stim_psth,pre_trial_psth,unitIds] = change_pre_stim_fr(file,500,500);
    
    pre_trial = mean(pre_trial_psth);
    pre_stim  = mean(pre_stim_psth);
    ratio{n}  = pre_stim./pre_trial;
    
end

figure,
for n = 1:length(ratio)
    plot(n,ratio{n},'ok'),hold on
end

title(subject)
xlim([0 length(ratio)+1])