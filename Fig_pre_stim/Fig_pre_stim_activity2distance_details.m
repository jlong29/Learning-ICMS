% the purspose of this code is to determine if the changes I've observed in the
% initial ICMS excitation across learning follows a spatial pattern with
% distance from the stimulating electrode. This also allows me to nicely
% visualize all the data for all animals

% Set window
win       = 1000;

% Subject IDs
subjects = {'R29','V1','V4','V7','V8'};
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig_pre_stim')
save('Fig_pre_stim_activity','win')

for q = 1:length(subjects)
    disp(subjects{q})
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
        [pd_base,pd_stim,unitIds] = change_pre_stim_fr(file,win);
        
        % Determine size of variable and allocate memory
        temp     = length(unitIds);

        % Data_matrix with the following columns:
        % 1) Difference in medians between paired samples
        % 2) Difference in variance between paired samples
        % 3) distance from stimulating electrode
        % 4) test(1) or control(0) hemisphere
        
        pre_stim2dist  = zeros(temp,4);
        
        for p = 1:temp
            
            % Calculate distance of unit from stimulating electrode
            loc = str2double(unitIds{p}(4:6));
            
            if loc <= 8
                pre_stim2dist(p,3) = .3 + .25*(loc-1);
                pre_stim2dist(p,4) = 1;
            elseif loc > 8 && loc <=16
                pre_stim2dist(p,3) = .3*sqrt(2) + .25*((loc-8)-1);
                pre_stim2dist(p,4) = 1;
            elseif loc < 25
                pre_stim2dist(p,3) = .3 + .25*(loc-16-1);
                pre_stim2dist(p,4) = 0;
            else
                pre_stim2dist(p,3) = .3*sqrt(2) + .25*((loc-16-8)-1);
                pre_stim2dist(p,4) = 0;
            end
            
            % Calculate difference ins paired sample medians (base2-base1 & stim
            % - base1)
            pre_stim2dist(p,1) = median(pd_stim(:,p)) - median(pd_base(:,p));
            pre_stim2dist(p,2) = var(pd_stim(:,p)) - var(pd_base(:,p));
            
        end

        eval([subjects{q} '_' date '_pre_stim2dist  = pre_stim2dist;'])

    end

    close all
    
    % Save it
    cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig_pre_stim')
    eval(['save(''Fig_pre_stim_activity'',''' subjects{q} '_*'',''-append'')'])
    eval(['clear(''' subjects{q} '*'')'])
    
end

system('%windir%\system32\rundll32.exe PowrProf.dll, SetSuspendState')

%% Now let's take a look at the relationship between ICMS response and distance
%%% from stimulating electrode for each subject
close all
clear

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig_pre_stim')
scrnsz = get(0,'screensize');

figure('Position',scrnsz)

% Subject IDs
subjects = {'R29*','V1*','V4*','V7*','V8*'};

for p = 1:length(subjects)
    load('Fig_pre_stim_activity',subjects{p})
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

            title([strrep(name,'_',' ') ': median(\Deltastim) - median(\Deltabase)']),
            ylabel('Difference of Medians')
            xlabel('Distance from stimulating electrode (mm)')
            axis square

        subplot(122),plot(sub(b,3),sub(b,2),'k'),hold on,
            plot(sub(b,3),sub(b,2),'ob','markerfacecolor','b'),
            ylim([ymin(2) ymax(2)])

            title([strrep(name,'_',' ') ': variance(\Deltastim) - variance(\Deltabase)'])
            ylabel('Difference of Variances')
            xlabel('Distance from stimulating electrode (mm)')
            axis square

            pause
            clf

    end
    
    eval(['clear(''' subjects{p} ''')'])
    
end

%% Example file for first two panels
% Let's grab the file data
file  = 'C:\Data\Plexon_Data\Test_Data\Learning\V8\V8_040109_spikes&Events1_16.mat';

% Look at paired differences between start trial samples and start
% trial/pre-stim sample
[pd_base,pd_stim,unitIds] = change_pre_stim_fr(file,1000);

upper = max([pd_stim(:);pd_base(:)]);
lower = min([pd_stim(:);pd_base(:)]);
edge = lower:upper;
