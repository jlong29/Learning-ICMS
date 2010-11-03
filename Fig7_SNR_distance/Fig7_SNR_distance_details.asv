% The purpose of this code is to determine whether or not there is a
% relationship between the magnitude of the SNR metrics and the distance of the
% recorded unit from the stimulating electrode.
% basic parameters
bin       = 10;
win       = 500;
ind       = win/bin;

% Subject IDs
subjects = {'R29','V1','V4','V7','V8'};
for q = 1:length(subjects)
    % Generate name of relevant subject directory
    folder = ['C:\Data\Plexon_Data\Test_Data\' subjects{q} '\' ...
            subjects{q} ' Spikes & Events'];
        
    cd(folder)
    
    % Grab the relevant files
    [file_list] = find_files(subjects{q},'spikes&Events.');

    % Generate SNR signal and internal noise models
    for n = 1:length(file_list)
        % Let's grab the file date
        file  = file_list{n};
        [temp,name] = fileparts(file);
        parse = findstr(name,'_');
        date  = name(parse(1)+1:parse(2)-1);

        % Generate PSTHs for each unit reference to ICMS
        [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file,'Event015',bin,win,win);

        %%%%%%%%%%%%%%%%%%%
        %%% SNR Metrics %%%
        %%%%%%%%%%%%%%%%%%%
        data     = all_psth;

        % Determine size of variable and allocate memory
        temp     = size(data,2);

        % Data_matrix with the following columns:
        % 1) signal model
        % 2) internal noise model
        % 3) distance from stimulating electrode
        % 4) test(1) or control(0) hemisphere

        snr2dist  = zeros(temp,4);

        for p = 1:temp

            % Calculate distance of unit from stimulating electrode
            loc = str2double(unitIds{p}(4:6));

            if loc <= 8
                snr2dist(p,3) = .3 + .25*(loc-1);
                snr2dist(p,4) = 1;
            elseif loc > 8 && loc <=16
                snr2dist(p,3) = .3*sqrt(2) + .25*((loc-8)-1);
                snr2dist(p,4) = 1;
            elseif loc < 25
                snr2dist(p,3) = .3 + .25*(loc-16-1);
                snr2dist(p,4) = 0;
            else
                snr2dist(p,3) = .3*sqrt(2) + .25*((loc-16-8)-1);
                snr2dist(p,4) = 0;
            end

            pre        = data(1:ind,p);
            post       = data(ind+1:end,p);

            % Calculate increase in signal relative to background (z-score)
            snr2dist(p,1) = (max(post) - mean(pre))/std(pre);

            % Calculate internal noise metric 
            meanb = mean(pre);

            % Generate baseline bin values below baseline average firing rate
            base  = meanb - pre;
            base  = base(base > 0);

            % Generate test bins values below baseline average firing rate
            test  = meanb - post;
            test  = test(test > 0);

            % Calculate increase in signal relative to background (t-test, unequal
            % samples and unequal variances)
            snr2dist(p,2) = (mean(test)-mean(base))/...
                          (sqrt(var(test)/length(test) + var(base)/length(base)));

        end

        eval([subjects{q} '_' date '_snr2dist  = snr2dist;'])

    end

    close all
    
    % Save it
    cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig7_SNR_distance')
    eval(['save(''Fig7_SNR_distance_test'',''' subjects{q} '_*'',''-append'')'])
    eval(['clear(''' subjects{q} '*'')'])
    
end

%% Now let's take a look at the relationship between ICMS response and distance
%%% from stimulating electrode for each subject
close all
clear

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig7_SNR_distance')

scrnsz = get(0,'screensize');
figure('Position',scrnsz)

% Subject IDs
subjects = {'R29*','V1*','V4*','V7*','V8*'};

for p = 1:length(subjects)
    load('Fig7_SNR_distance_test',subjects{p})
    vars = who(subjects{p});
    
    % Grab max and min for signal and internal noise
    maxs  = zeros(length(vars),2);
    mins  = zeros(length(vars),2);

    for n = 1:length(vars)
        temp    = eval(vars{n});
        temp(isinf(temp)) = NaN;

        maxs(n,:) = nanmax(temp(:,1:2));
        mins(n,:) = nanmin(temp(:,1:2));
    end

    % Calculate the max and min for each metric
    ymin = floor(min(mins));
    ymax = ceil(max(maxs));

    for n = 1:length(vars)

        name = vars{n};
        temp = eval(vars{n});
        
        if strcmp(subjects{p},'R29*')
            ind  = ~logical(temp(:,4));
        else
            ind  = logical(temp(:,4));
        end
        
        sub = temp(ind,:);
        [a,b]=sort(sub(:,3),1);

        subplot(121),plot(sub(b,3),sub(b,1),'k'),hold on,
            plot(sub(b,3),sub(b,1),'ob','markerfacecolor','b'),
            ylim([ymin(1) ymax(1)])

            title([strrep(name,'_',' ') ': post-stim excitation']),
            ylabel('z-score')
            xlabel('Distance from stimulating electrode (mm)')
            axis square

        subplot(122),plot(sub(b,3),sub(b,2),'k'),hold on,
            plot(sub(b,3),sub(b,2),'ob','markerfacecolor','b'),
            ylim([ymin(2) ymax(2)])

            title([strrep(name,'_',' ') ': post-stim inhibition'])
            ylabel('t-value')
            xlabel('Distance from stimulating electrode (mm)')
            axis square

            pause
            clf

    end
    
    eval(['clear(''' subjects{p} ''')'])
    
end

%% Let's see if there is a trend when data is grouped across all subjects around
%%% the entrance into stage 2.
close all
clear

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig7_SNR_distance')

scrnsz = get(0,'screensize');
figure('Position',scrnsz)

% Subject IDs
subjects = {'R29*','V1*','V4*','V7*','V8*'};
load('Fig7_SNR_distance_test_sub2')

% Grab the data
vars1 = who(subjects{1});
vars2 = who(subjects{2});
vars3 = who(subjects{3});
vars4 = who(subjects{4});
vars5 = who(subjects{5});
    
% Grab max and min for signal and internal noise
maxs  = zeros(length(vars1),2);
mins  = zeros(length(vars1),2);

for n = 1:length(vars1)
    temp1    = eval(vars1{n});
    temp2    = eval(vars2{n});
    temp3    = eval(vars3{n});
    temp4    = eval(vars4{n});
    temp5    = eval(vars5{n});
    
    temp1(isinf(temp1)) = NaN;
    temp2(isinf(temp2)) = NaN;
    temp3(isinf(temp3)) = NaN;
    temp4(isinf(temp4)) = NaN;
    temp5(isinf(temp5)) = NaN;
    
    maxs(n,:) = nanmax([temp1(:,1:2);temp2(:,1:2);temp3(:,1:2);temp4(:,1:2);temp5(:,1:2)]);
    mins(n,:) = nanmin([temp1(:,1:2);temp2(:,1:2);temp3(:,1:2);temp4(:,1:2);temp5(:,1:2)]);
end

% Calculate the max and min for each metric
ymin = floor(min(mins));
ymax = ceil(max(maxs));

% temp to correct for outlier
ymax(1) = 50;

for n = 1:length(vars1)
    
    % Data for each subject
    temp1 = eval(vars1{n});
    temp2 = eval(vars2{n});
    temp3 = eval(vars3{n});
    temp4 = eval(vars4{n});
    temp5 = eval(vars5{n});
    
    ind1  = ~logical(temp1(:,4));
    ind2  = logical(temp2(:,4));
    ind3  = logical(temp3(:,4));
    ind4  = logical(temp4(:,4));
    ind5  = logical(temp5(:,4));
    
    % Test Hemisphere data for each subject
    sub1  = temp1(ind1,:);
    [a,b] = sort(sub1(:,3),1);
    sub1  = sub1(b,:);
    
    sub2  = temp2(ind2,:);
    [a,b] = sort(sub2(:,3),1);
    sub2  = sub2(b,:);
    
    sub3  = temp3(ind3,:);
    [a,b] = sort(sub3(:,3),1);
    sub3  = sub3(b,:);
    
    sub4  = temp4(ind4,:);
    [a,b] = sort(sub4(:,3),1);
    sub4  = sub4(b,:);
    
    sub5  = temp5(ind5,:);
    [a,b] = sort(sub5(:,3),1);
    sub5  = sub5(b,:);
    
    subplot(2,3,n),plot(sub1(:,3),sub1(:,1),'+k'),hold on,
        plot(sub2(:,3),sub2(:,1),'ok'),
        plot(sub3(:,3),sub3(:,1),'sk'),
        plot(sub4(:,3),sub4(:,1),'vk'),
        plot(sub5(:,3),sub5(:,1),'*k'),
        
        % fit a line to the group data
        h = get(gca,'children');
        x = cell2mat(get(h,'xdata')');
        y = cell2mat(get(h,'ydata')');
        p = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        x = 0:.2:2.25;
        y = polyval(p,x);
        plot(x,y,'--k')
        
        ylim([ymin(1) ymax(1)])
        xlim([0 2.25])
        ylabel('z-score')
        xlabel('Distance from stimulating electrode (mm)','fontsize',8)
        axis square
        
        if n == 1
            title('ICMS Signal: Stage1: Confused','fontsize',8)
            legend('S1','S2','S3','S4','S5','location','northwest')
        elseif n == 2
            title('ICMS Signal: Stage2: Learning','fontsize',8)
        else
            title('ICMS Signal: Stage3: Acquired','fontsize',8)
        end
        
    subplot(2,3,n+3),plot(sub1(:,3),sub1(:,2),'+k'),hold on,
        plot(sub2(:,3),sub2(:,2),'ok'),
        plot(sub3(:,3),sub3(:,2),'sk'),
        plot(sub4(:,3),sub4(:,2),'vk'),
        plot(sub5(:,3),sub5(:,2),'*k'),
        
        % fit a line to the group data
        h = get(gca,'children');
        x = cell2mat(get(h,'xdata')');
        y = cell2mat(get(h,'ydata')');
        p = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        x = 0:.2:2.25;
        y = polyval(p,x);
        plot(x,y,'--k')
        
        ylim([ymin(2) ymax(2)])
        xlim([0 2.25])
        title(['Day ' num2str(n) ': internal noise'])
        ylabel('t-value')
        xlabel('Distance from stimulating electrode (mm)','fontsize',8)
        axis square

        if n == 1
            title('ICMS Internal Noise: Stage1: Confused','fontsize',8)
            legend('S1','S2','S3','S4','S5','location','northwest')
        elseif n == 2
            title('ICMS Internal Noise: Stage2: Learning','fontsize',8)
        else
            title('ICMS Internal Noise: Stage3: Acquired','fontsize',8)
        end
        
end

set(gcf,'color','w')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 7.6])

