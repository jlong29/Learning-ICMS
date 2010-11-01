%%%%%%%%%%%%%%%%%%%%
%%% Example Units %%%
%%%%%%%%%%%%%%%%%%%%%
%% Grab Example units from .plx file
save_folder = 'C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures';

% R29
cd('E:\John_backup\John Plexon Data\r29\')

units  = cell(1,3);
units1 = {'20a';'25a';'26b'};
units2 = {'17a';'17b';'19a';'22a';'30a'};
units3 = {'17a';'18b';'28a'};

units{1} = units1;
units{2} = units2;
units{3} = units3;

file_names = {'R29_050808_rec1HD_blockICMS_clean.plx';
              'R29_051008_rec1HD_blockICMS_clean.plx';
              'R29_051208_rec1HD_blockICMS_clean.plx'};
tic
for p = 1:length(units)
    p
    for q = 1:length(units{p})
        q
        if strcmpi(units{p}{q}(3),'a')
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 1);
            eval(['S' num2str(p) '_R29_sig0' units{p}{q} ' = wave;'])
        else
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 2);
            eval(['S' num2str(p) '_R29_sig0' units{p}{q} ' = wave;'])
        end
    end
    save([save_folder '\' file_names{p}(1:end-4) '.mat'],'-regexp',['S' num2str(p)])
    clear S*
end
toc

% V7
cd('E:\John_backup\John Plexon Data\V7\')

units  = cell(1,3);
units1 = {'02a';'11a'};
units2 = {'02a';'04a';'07a';'08a';'16a'};
units3 = {'13a'};

units{1} = units1;
units{2} = units2;
units{3} = units3;

file_names = {'V7_032709_ICMS_clean.plx';
              'V7_032909_baseline_and_ICMS_clean.plx';
              'V7_033109_ICMS_clean.plx'};
tic
for p = 1:length(units)
    p
    for q = 1:length(units{p})
        q
        if strcmpi(units{p}{q}(3),'a')
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 1);
            eval(['S' num2str(p) '_V7_sig0' units{p}{q} ' = wave;'])
        else
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 2);
            eval(['S' num2str(p) '_V7_sig0' units{p}{q} ' = wave;'])
        end
    end
    save([save_folder '\' file_names{p}(1:end-4) '.mat'],'-regexp',['S' num2str(p)])
    clear S*
end
toc

% V8
cd('E:\John_backup\John Plexon Data\V8\')

units  = cell(1,3);
units1 = {'11a'};
units2 = {'05a';'08a';'11a';'13a'};
units3 = {'04a';'11a'};

units{1} = units1;
units{2} = units2;
units{3} = units3;

file_names = {'V8_032909_ICMS_clean.plx';
              'V8_033109_baseline_and_ICMS_clean.plx';
              'V8_040109_ICMS_clean.plx'};
tic
for p = 1:length(units)
    p
    for q = 1:length(units{p})
        q
        if strcmpi(units{p}{q}(3),'a')
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 1);
            eval(['S' num2str(p) '_V8_sig0' units{p}{q} ' = wave;'])
        else
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 2);
            eval(['S' num2str(p) '_V8_sig0' units{p}{q} ' = wave;'])
        end
    end
    save([save_folder '\' file_names{p}(1:end-4) '.mat'],'-regexp',['S' num2str(p)])
    clear S*
end
toc

% V1
cd('E:\John_backup\John Plexon Data\V1\')

units  = cell(1,3);
units1 = {'03a';'16a'};
units2 = {'01a';'16a'};
units3 = {'02a';'05a'};

units{1} = units1;
units{2} = units2;
units{3} = units3;

file_names = {'V1_092408_1HD_ICMS_tr_clean.plx';
              'V1_092608_ICMS_clean.plx';
              'V1_093008_ICMS_clean.plx'};
tic
for p = 1:length(units)
    p
    for q = 1:length(units{p})
        q
        if strcmpi(units{p}{q}(3),'a')
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 1);
            eval(['S' num2str(p) '_V1_sig0' units{p}{q} ' = wave;'])
        else
            % 1-based channel number system. Units a-d are 1-4
            [n, npw, ts, wave] = plx_waves_v(file_names{p}, str2double(units{p}{q}(1:2)), 2);
            eval(['S' num2str(p) '_V1_sig0' units{p}{q} ' = wave;'])
        end
    end
    save([save_folder '\' file_names{p}(1:end-4) '.mat'],'-regexp',['S' num2str(p)])
    clear S*
end
toc

%% Divide Waveforms into beginning and end of session and clean
% I use an L2 norm criteria for clean
% Go to the example cell folder
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')

R29_files = dir('R29*');
V7_files  = dir('V7*');
V8_files  = dir('V8*');
V1_files  = dir('V1*');

% Organize files by stage
files_S1 = {R29_files(1).name; V7_files(1).name; V8_files(1).name; V1_files(1).name}; 
files_S2 = {R29_files(2).name; V7_files(2).name; V8_files(2).name; V1_files(2).name}; 
files_S3 = {R29_files(3).name; V7_files(3).name; V8_files(3).name; V1_files(3).name}; 

% For each stage find the number of units (for figure layout)
for n = 1:length(files_S1)
    load(files_S1{n})
end

numS1 = length(who('S1*'));
clear S1*

for n = 1:length(files_S2)
    load(files_S2{n})
end

numS2 = length(who('S2*'));
clear S2*

for n = 1:length(files_S3)
    load(files_S3{n})
end

numS3 = length(who('S3*'));
clear S3*

%% Stage 1 Example Units
scrnsz = get(0,'screensize');
figure('Position',scrnsz,'Name','Stage 1 of Learning ICMS, Example Units')

% Generate figure of example units during S1
fig_ind = 1;

for n = 1:length(files_S1)
    load(files_S1{n})
    unitIds = who('S1*');
    
    for p = 1:length(unitIds)
        
        eval(['start = ' unitIds{p} '(1:1000,:);'])
        eval(['fin   = ' unitIds{p} '(end-999:end,:);'])
                
        meanst     = mean(start);
        strss      = sum((start - repmat(meanst,1000,1)).^2,2);
        
        meanfin    = mean(fin);
        finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);
        
        [c,stind]   = sort(strss);
        [d,finind] = sort(finrss);
        h = subplot(2,numS1,fig_ind); plot(start(stind(1:200),:)','k')
            xlim([1 32])
            temp  = logical(unitIds{p} == '_');
            label = unitIds{p};
            label(temp) = ' ';
            title(label)
            y      = get(h,'ylim');
            set(h,'ylim',y)
            ticks  = get(h,'ytick');
            set(h,'ytick',[ticks(1) 0 ticks(end)])
            ticks  = get(h,'ytick');
            lticks = get(h,'yticklabel');
            axis square
        
        j = subplot(2,numS1,fig_ind+numS1); plot(fin(finind(1:200),:)','k')
            xlim([1 32])
            set(j,'ylim',y)
            set(j,'ytick',ticks)
            set(j,'yticklabel',lticks)
            axis square
            
            drawnow
            
        fig_ind    = fig_ind +1;
    end
    clear S1*
end

%% Stage 2 Example Units
scrnsz = get(0,'screensize');
figure('Position',scrnsz,'Name','Stage 2 of Learning ICMS, Example Units')

% Generate figure of example units during S1
fig_ind = 1;

for n = 1:length(files_S2)
    load(files_S2{n})
    unitIds = who('S2*');
    
    for p = 1:length(unitIds)
        
        eval(['start = ' unitIds{p} '(1:1000,:);'])
        eval(['fin   = ' unitIds{p} '(end-999:end,:);'])
                
        meanst     = mean(start);
        strss      = sum((start - repmat(meanst,1000,1)).^2,2);
        
        meanfin    = mean(fin);
        finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);
        
        [c,stind]   = sort(strss);
        [d,finind] = sort(finrss);
        h = subplot(4,8,fig_ind); plot(start(stind(1:200),:)','k')
            xlim([1 32])
            temp  = logical(unitIds{p} == '_');
            label = unitIds{p};
            label(temp) = ' ';
            title(label)
            y      = get(h,'ylim');
            set(h,'ylim',y)
            ticks  = get(h,'ytick');
            set(h,'ytick',[ticks(1) 0 ticks(end)])
            ticks  = get(h,'ytick');
            lticks = get(h,'yticklabel');
            axis square
            
        j = subplot(4,8,fig_ind+numS2/2); plot(fin(finind(1:200),:)','k')
            xlim([1 32])
            set(j,'ylim',y)
            set(j,'ytick',ticks)
            set(j,'yticklabel',lticks)
            axis square
            
            drawnow
            
        if fig_ind == 8
            fig_ind =17;
        else
            fig_ind    = fig_ind +1;
        end
    end
    clear S2*
end

%% Stage 3 Example Units
scrnsz = get(0,'screensize');
figure('Position',scrnsz,'Name','Stage 3 of Learning ICMS, Example Units')

% Generate figure of example units during S1
fig_ind = 1;

for n = 1:length(files_S3)
    load(files_S3{n})
    unitIds = who('S3*');
    
    for p = 1:length(unitIds)
        
        eval(['start = ' unitIds{p} '(1:1000,:);'])
        eval(['fin   = ' unitIds{p} '(end-999:end,:);'])
                
        meanst     = mean(start);
        strss      = sum((start - repmat(meanst,1000,1)).^2,2);
        
        meanfin    = mean(fin);
        finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);
        
        [c,stind]   = sort(strss);
        [d,finind] = sort(finrss);
        h = subplot(2,numS3,fig_ind); plot(start(stind(1:200),:)','k')
            xlim([1 32])
            temp  = logical(unitIds{p} == '_');
            label = unitIds{p};
            label(temp) = ' ';
            title(label)
            y      = get(h,'ylim');
            set(h,'ylim',y)
            ticks  = get(h,'ytick');
            set(h,'ytick',[ticks(1) 0 ticks(end)])
            ticks  = get(h,'ytick');
            lticks = get(h,'yticklabel');
            axis square
            
        j = subplot(2,numS3,fig_ind+numS3); plot(fin(finind(1:200),:)','k')
            xlim([1 32])
            set(j,'ylim',y)
            set(j,'ytick',ticks)
            set(j,'yticklabel',lticks)
            axis square
            
            drawnow
        
        fig_ind = fig_ind +1;
        
    end
    clear S3*
end

%% Selecting Among the finalists
clear
scrnsz = get(0,'screensize');
figure('Position',scrnsz)

% I use an L2 norm criteria for clean
% Go to the example cell folder
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')

R29_files = dir('R29*');
V7_files  = dir('V7*');
V8_files  = dir('V8*');
V1_files  = dir('V1*');

% Organize files by stage
files_S1 = {R29_files(1).name; V7_files(1).name; V8_files(1).name; V1_files(1).name}; 
files_S2 = {R29_files(2).name; V7_files(2).name; V8_files(2).name; V1_files(2).name}; 
files_S3 = {R29_files(3).name; V7_files(3).name; V8_files(3).name; V1_files(3).name}; 

% For each stage find the number of units (for figure layout)
for n = 1:length(files_S1)
    load(files_S1{n})
end

numS1 = length(who('S1*'));
clear S1*

for n = 1:length(files_S2)
    load(files_S2{n})
end

numS2 = length(who('S2*'));
clear S2*

for n = 1:length(files_S3)
    load(files_S3{n})
end

numS3 = length(who('S3*'));
clear S3*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's take a closer look at Stage 2 units %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject_dirs    = {'C:\Data\Plexon_Data\Test_Data\Learning\R29'; 'C:\Data\Plexon_Data\Test_Data\Learning\V7';...
    'C:\Data\Plexon_Data\Test_Data\Learning\V8';'C:\Data\Plexon_Data\Test_Data\Learning\V1'};
source_files_S2 = {'R29_051008_spike&Events17_32.mat'; 'V7_032909_spikes&Events1_16.mat'; ...
    'V8_033109_spikes&Events1_16.mat'; 'V1_092608_spikes&Events_1_16.mat'};

for i = 1:length(files_S1)
    
    cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig2_Example_Units')
    load(files_S2{i})
    wf_list = who('S*');
    
    for j = 1:length(wf_list)
        
        unit    = wf_list{j};

        cd(subject_dirs{i})
        load(source_files_S2{i},unit(end-6:end))
        sp_list = who('sig*');
        
        eval(['spikes = ' sp_list{j} ';'])
        eval(['wf     = ' wf_list{j} ';'])

        % Basic parameters
        %%%%%%%%%%%%%%%%%%
        samples = -200:32:800;  % us time axis for waveforms
        edge    = 0:.001:.1;    % binning for ISI calculation
        win     = 800;          % window for calculating firing rate in seconds
        bin     = 10;           % bin for above window in seconds

        start   = wf(1:1000,:);
        isi1    = diff(spikes(1:5001));

        fin     = wf(end-999:end,:);
        isi2    = diff(spikes(end-5000:end,:));

        meanst     = mean(start);
        strss      = sum((start - repmat(meanst,1000,1)).^2,2);

        meanfin    = mean(fin);
        finrss     = sum((fin - repmat(meanfin,1000,1)).^2,2);

        [c,stind]  = sort(strss);
        [d,finind] = sort(finrss);

        % Look at early session waveforms
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         h = subplot(1,4,1); plot(start(stind(1:200),:)','k')
        h = subplot(1,4,1); plot(start(1:200,:)','k')
            xlim([1 32])
            temp  = logical(unit == '_');
            unit(temp) = ' ';
            title(unit)
            y      = get(h,'ylim');
            set(h,'ylim',y)
            ticks  = get(h,'ytick');
            set(h,'ytick',[ticks(1) 0 ticks(end)])
            ticks  = get(h,'ytick');
            lticks = get(h,'yticklabel');
            axis square

        % Look at late session waveforms
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         j = subplot(1,4,2); plot(fin(finind(1:200),:)','k')
        j = subplot(1,4,2); plot(fin(1:200,:)','k')
            xlim([1 32])
            set(j,'ylim',y)
            set(j,'ytick',ticks)
            set(j,'yticklabel',lticks)
            axis square

        % Look at early and late session ISIs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [h,p] = kstest2(isi1,isi2,.01,'unequal');
        subplot(1,4,3),plot(edge,histc(isi1,edge),'color',[0.5412 0.5490 0.4510],'linewidth',2),hold on,
                       plot(edge,histc(isi2,edge),'k','linewidth',2)
                       xlim([-.005 edge(end)])
                       title(['KS p = ' sprintf('%3.2f',p)])
                       axis square

        % Look at early and late session firing rate of unit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        time_line1 = spikes(1):bin:spikes(1) + win;

        unit_fr1   = zeros(1,length(time_line1)-1);
        data1      = spikes;

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
        % time_line2 = 2200:bin:2560;
        time_line2 = spikes(end) - win:bin:spikes(end);

        unit_fr2   = zeros(1,length(time_line2)-1);
        data1      = spikes;

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

        subplot(1,4,4), bar(1,mean1,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
                        bar(2,mean2,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),
                        errorbar([1 2],[mean1 mean2],[-std1 -std2],[std1 std2],'color','k','linestyle','none','linewidth',3)
                        set(gca,'xtick',[1 2])
                        set(gca,'xticklabel',{'Early';'Late'})
                        title([{['Mean1 = ' num2str(mean1)]};{['Variance1 = ' num2str(var(unit_fr1))]};
                               {['Mean2 = ' num2str(mean2)]};{['Variance2 = ' num2str(var(unit_fr2))]}])
                        axis square
        pause
        clf
    end
    keep files_S2 subject_dirs source_files_S2
end
