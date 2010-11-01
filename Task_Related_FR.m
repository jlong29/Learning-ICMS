% Comparing Unit firing rates during baseline, tone task, and ICMS task
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Global Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%
win        = 20000; % samples in milliseconds

%%%%%%%%%%%%%%%
%%% Subject %%%
%%%%%%%%%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\V1\V1 Spikes & Events')

% Grab all files with reference to ICMS task
ICMS_files = find_files('V1','spikes&Events.');

% Only include relevant files (arbitrary by subject)
ICMS_files = ICMS_files(1:end-1);

% Let's look for separate '_tone' files that contain tone and baseline task
AUX_files  = cell(size(ICMS_files));

for n = 1:length(ICMS_files)
    [path,file]   = fileparts(ICMS_files{n});
    % Get rid of full path name for ICMS files
    ICMS_files{n} = file;
    
    chk           = find_files(file,'_tone');
    
    if ~isempty(chk)
        AUX_files{n} = chk{:};
    else
        AUX_files{n} = '';
    end
    
end

% In the case of corresponding '_tone' files, sometimes these contain just 
% the baseline and sometimes these contain the baseline and tone task data
% 1) If '_tone' file has no Event011, then just baseline
% 2) If no '_tone' file and Event011(1) > 300, then file contains baseline
%       and tone task data
% 3) If no '_tone' file and post 'Behavior_preproc' only Event011test or
%       Event011train has > 5 events, then no tone task data.

% Let's find out whether the baseline and Tone data exist and where it is
% located

Tone_files = cell(ICMS_files);
Base_files = cell(ICMS_files);

for n = 1:length(ICMS_files)
    
    % Examine '_tone' file if it exists
    if ~isempty(AUX_files{n})
        [path,file] = fileparts(AUX_files{n});
        
        % Check for Event011
        chk         = who('-file', file,'-regexp','Event011');
        if isempty(chk)
            % '_tone' file is just baseline data
            Base_files{n} = file;
            % Tone data might be in ICMS_file
            Tone_files{n} = ICMS_files{n};
        else
            % '_tone' file is baseline and tone data
            Base_files{n} = file;
            Tone_files{n} = file;
        end
        
    % Otherwise, check for baseline and tone data is main ICMS file
    else
        chk         = load(ICMS_files{n},'Event011');
        % Confirm baseline data
        if chk.Event011(1) > 300
            Base_files{n} = ICMS_files{n};
        else
            Base_files{n} = '';
        end
        
        % Set ICMS_file as candidate for tone data
        Tone_files{n} = ICMS_files{n};
    end
end

% Having Separated all the files, now it is time to compare firing rates
% between task conditions
for q = 1:length(ICMS_files)
    q
    chk1 = strcmp(ICMS_files{q},Tone_files{q});
    chk2 = strcmp(ICMS_files{q},Base_files{q});
    chk3 = strcmp(Tone_files{q},Base_files{q});
    
    % Case 1:
    % if chk1 =0 and chk2 = 0 and chk3 = 1, then baseline and tone data
    % are in a separate file
    if chk1 == 0 && chk2 == 0 && chk3 == 1
        %%% Grab Baseline and Tone data %%%
        load(Tone_files{q})
        
        % Generate formatted data
        data_matrix
        
        % Save associated unit IDs
        base_list = sp_list;
        tone_list = sp_list;
        
        % partition between baseline and tone data
        part = find(time_line < Event011(1),1,'last');
        
        % Separate Data accordingly
        base_data = spikes(1:part-1,:);
        tone_data = spikes(part:end,:);
        
        time_line = time_line(part:end);
        
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        % check how trials were logged
        if length(event.Event011train) > 5
            class = 'train';
        else
            class = 'test';
        end
        
        [tone_data] = inter_trial_data(tone_data,time_line,event,class);
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Grab ICMS data %%%
        load(ICMS_files{q})
        
        % Generate formatted data
        data_matrix
        
        % Save associated unit IDs and data
        icms_list = sp_list;
        icms_data = spikes;
        
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        
        [icms_data] = inter_trial_data(icms_data,time_line,event,'test');
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Only look at units that are common among the files %%%
        %%% Compare equal amounts of data                      %%%
        [ind1,ind2] = list_comp(icms_list,base_list);
        [ind1,ind3] = list_comp(icms_list(ind1),tone_list);
        
        [val] = min([size(icms_data,1) size(base_data,1) size(tone_data,1)]);
        
        % Update data
        base_list   = base_list(ind2);
        base_data   = base_data(1:val,ind2);
        
        tone_list   = tone_list(ind3);
        tone_data   = tone_data(1:val,ind3);
        
        icms_list   = icms_list(ind1);
        icms_data   = icms_data(1:val,ind1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Code for calculating firing rate after all data corrections
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        resid    = mod(val,win);
        
        % Determine number of non-overlapping segments of length win
        N        = (val-resid)/win;
        
        % Set memory for ensemble firing rates
        base_fr  = zeros(N, sum(ind1));
        tone_fr  = zeros(N, sum(ind1));
        icms_fr  = zeros(N, sum(ind1));
        
        for p    = 1:N
            base_fr(p,:) = sum(base_data(win*(p-1)+1:win*p,:));
            tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
            icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
        end
        
        % Convert to spike per second
        base_fr   = [mean(base_fr.*(1000/win))' std(base_fr.*(1000/win))'];
        tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
        icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
    % Case 2:
    % if chk1 = 0 and chk2 = 0 and chk3 == 0, tone and icms data are in
    % different files but there is no baseline data
    elseif chk1 == 0 && chk2 == 0 && chk3 == 0
        %%% Grab Baseline and Tone data %%%
        load(Tone_files{q})
        
        % Generate formatted data
        data_matrix
        
        % Save associated unit IDs
        base_list = '';
        tone_list = sp_list;
        
        % Separate Data accordingly
        tone_data = spikes;
        
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        % check how trials were logged
        if length(event.Event011train) > 5
            class = 'train';
        else
            class = 'test';
        end
        
        [tone_data] = inter_trial_data(tone_data,time_line,event,class);
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Grab ICMS data %%%
        load(ICMS_files{q})
        
        % Generate formatted data
        data_matrix
        
        % Save associated unit IDs and data
        icms_list = sp_list;
        icms_data = spikes;
        
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        
        [icms_data] = inter_trial_data(icms_data,time_line,event,'test');
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Only look at units that are common among the files %%%
        %%% Compare equal amounts of data                      %%%
        [ind1,ind2] = list_comp(icms_list,tone_list);
        
        [val] = min([size(icms_data,1) size(tone_data,1)]);
        
        % Update data
        tone_list   = tone_list(ind2);
        tone_data   = tone_data(1:val,ind2);
        
        icms_list   = icms_list(ind1);
        icms_data   = icms_data(1:val,ind1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Code for calculating firing rate after all data corrections
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        resid    = mod(val,win);
        
        % Determine number of non-overlapping segments of length win
        N        = (val-resid)/win;
        
        % Set memory for ensemble firing rates
        tone_fr  = zeros(N, sum(ind1));
        icms_fr  = zeros(N, sum(ind1));
        
        for p    = 1:N
            tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
            icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
        end
        
        
        % Convert to spike per second
        tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
        icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Case 3:
    % if chk1 = 1 and chk2 = 0 and chk3 == 0, tone and icms data are
    % perhaps in the same file and baseline data may or may not exist in a
    % separate file
    elseif chk1 == 1 && chk2 == 0 && chk3 == 0
        %%% Grab Baseline data, if it exists %%%
        if isempty(Base_files{q})
            base_list = '';
            base_data = '';
        else
            load(Base_files{q})
            data_matrix
            
            base_list = sp_list;
            base_data = spikes;
            
            % A bit of Clean Up
            clear event Eve* spikes time_line sp_list
        end
        
        %%% Grab ICMS and TONE data %%%
        load(ICMS_files{q})
        
        % Generate formatted data
        data_matrix
        
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        
        % Grab ICMS data
        icms_list             = sp_list;
        [icms_data] = inter_trial_data(spikes,time_line,event,'test');
            
        % Grab Tone data if it exists
        if length(event.Event011train) > 5
            tone_list             = sp_list;
            [tone_data] = inter_trial_data(spikes,time_line,event,'train');
        else
            tone_list             = '';
            tone_data             = '';
        end
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Only look at units that are common among the files %%%
        %%% Compare equal amounts of data                      %%%
        if isempty(base_list) && ~isempty(tone_list)
            [ind1,ind2] = list_comp(icms_list,tone_list);
            
            [val] = min([size(icms_data,1) size(tone_data,1)]);
            
            % Update data
            tone_list   = tone_list(ind2);
            tone_data   = tone_data(1:val,ind2);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            tone_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif isempty(base_list) && isempty(tone_list)
            continue
        elseif ~isempty(base_list) && isempty(tone_list)
            [ind1,ind2] = list_comp(icms_list,base_list);
            
            [val] = min([size(icms_data,1) size(base_data,1)]);
            
            % Update data
            base_list   = base_list(ind2);
            base_data   = base_data(1:val,ind2);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            base_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                base_fr(p,:) = sum(base_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            base_fr   = [mean(base_fr.*(1000/win))' std(base_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            [ind1,ind2] = list_comp(icms_list,base_list);
            [ind1,ind3] = list_comp(icms_list(ind1),tone_list);
            
            [val] = min([size(icms_data,1) size(base_data,1) size(tone_data,1)]);
            
            % Update data
            base_list   = base_list(ind2);
            base_data   = base_data(1:val,ind2);
            
            tone_list   = tone_list(ind3);
            tone_data   = tone_data(1:val,ind3);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            base_fr  = zeros(N, sum(ind1));
            tone_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                base_fr(p,:) = sum(base_data(win*(p-1)+1:win*p,:));
                tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            base_fr   = [mean(base_fr.*(1000/win))' std(base_fr.*(1000/win))'];
            tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    % Case 4:
    % if chk1 = 1 and chk2 = 1 and chk3 == 1, everything might be in one
    % file
    elseif chk1 == 1 && chk2 == 1 && chk3 == 1
        
        %%% Load main file %%%
        load(ICMS_files{q})
        
        % Generate formatted data
        data_matrix
        
        %%% Grab Baseline data, if it exists %%%
        if Event011 < 300
            base_list = '';
            base_data = '';
        else
            % partition between baseline and other data
            part = find(time_line < Event011(1),1,'last');
            
            base_list = sp_list;
            base_data = spikes(1:part,:);
        end
      
        %%% Only look at task related data, not cue related data %%%
        event   = structpack(who('Eve*'));
        [event] = Behavior_preproc(event);
        
        % Grab ICMS data
        icms_list             = sp_list;
        [icms_data] = inter_trial_data(spikes,time_line,event,'test');
            
        % Grab Tone data if it exists
        if length(event.Event011train) > 5
            tone_list             = sp_list;
            [tone_data] = inter_trial_data(spikes,time_line,event,'train');
        else
            tone_list             = '';
            tone_data             = '';
        end
        
        % A bit of Clean Up
        clear event Eve* spikes time_line sp_list
        
        %%% Only look at units that are common among the files %%%
        %%% Compare equal amounts of data                      %%%
        if isempty(base_list) && ~isempty(tone_list)
            [ind1,ind2] = list_comp(icms_list,tone_list);
            
            [val] = min([size(icms_data,1) size(tone_data,1)]);
            
            % Update data
            tone_list   = tone_list(ind2);
            tone_data   = tone_data(1:val,ind2);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            tone_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif isempty(base_list) && isempty(tone_list)
            continue
        elseif ~isempty(base_list) && isempty(tone_list)
            [ind1,ind2] = list_comp(icms_list,base_list);
            
            [val] = min([size(icms_data,1) size(base_data,1)]);
            
            % Update data
            base_list   = base_list(ind2);
            base_data   = base_data(1:val,ind2);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            base_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                base_fr(p,:) = sum(base_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            base_fr   = [mean(base_fr.*(1000/win))' std(base_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            [ind1,ind2] = list_comp(icms_list,base_list);
            [ind1,ind3] = list_comp(icms_list(ind1),tone_list);
            
            [val] = min([size(icms_data,1) size(base_data,1) size(tone_data,1)]);
            
            % Update data
            base_list   = base_list(ind2);
            base_data   = base_data(1:val,ind2);
            
            tone_list   = tone_list(ind3);
            tone_data   = tone_data(1:val,ind3);
            
            icms_list   = icms_list(ind1);
            icms_data   = icms_data(1:val,ind1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Code for calculating firing rate after all data corrections
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            resid    = mod(val,win);
            
            % Determine number of non-overlapping segments of length win
            N        = (val-resid)/win;
            
            % Set memory for ensemble firing rates
            base_fr  = zeros(N, sum(ind1));
            tone_fr  = zeros(N, sum(ind1));
            icms_fr  = zeros(N, sum(ind1));
            
            for p    = 1:N
                base_fr(p,:) = sum(base_data(win*(p-1)+1:win*p,:));
                tone_fr(p,:) = sum(tone_data(win*(p-1)+1:win*p,:));
                icms_fr(p,:) = sum(icms_data(win*(p-1)+1:win*p,:));
            end
            
            % Convert to spike per second
            base_fr   = [mean(base_fr.*(1000/win))' std(base_fr.*(1000/win))'];
            tone_fr   = [mean(tone_fr.*(1000/win))' std(tone_fr.*(1000/win))'];
            icms_fr   = [mean(icms_fr.*(1000/win))' std(icms_fr.*(1000/win))'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end
    
    %%%%%%%%%%%%%%
    %%% Figure %%%
    %%%%%%%%%%%%%%
    % Set dots per inch
    dpi    = 96;
    scrnsz = get(0,'ScreenSize');
    
    figure('units','pixels','Position',scrnsz)
    % From left to right
    H = [1.75 5.5 9.25 ];
    
    % From bottom to top
    V = 3.5;
    
    % Layout of axes
    pos      = [H(1)*dpi  V(1)*dpi  2.5*dpi 2.5*dpi;
                H(2)*dpi  V(1)*dpi  2.5*dpi 2.5*dpi;
                H(3)*dpi  V(1)*dpi  2.5*dpi 2.5*dpi];
    
    % Baseline v. Tone data plot
    h1 = axes('units','pixels','position',pos(1,:));
        if isempty(base_list) || isempty(tone_list)
            text(0.1,0.5,'Data not Available for Comparison')
            set(h1,'xtick',[],'xticklabel','','ytick',[],'yticklabel','')
            title('Baseline v. Tone Data Plot')
        else
            ul = ceil(max([base_fr(:,1)+tone_fr(:,2); tone_fr(:,1)+base_fr(:,2)]));
            set(h1,'xlim',[0 ul],'ylim',[0 ul])
            line([0 ul],[0 ul],'color','k','linewidth',2,'linestyle','--')
            line([base_fr(:,1)'; base_fr(:,1)'],[tone_fr(:,1)'-base_fr(:,2)'; tone_fr(:,1)'+ base_fr(:,2)'],'color','k','linewidth',2)
            line([base_fr(:,1)'-tone_fr(:,2)'; base_fr(:,1)'+ tone_fr(:,2)'],[tone_fr(:,1)'; tone_fr(:,1)'],'color','k','linewidth',2)
            hold on
            
            % Separate Units by channels 1-16 and 17-32
            ind = zeros(length(base_list),1);
            for j = 1:length(base_list)
                temp = str2double(base_list{j}(5:6));
                if temp < 17
                    ind(j) = 1;
                end
            end
            ind = logical(ind);
            
            % Now plot it out
            plot(base_fr(ind,1),tone_fr(ind,1),'ok','markerfacecolor','c'),
            plot(base_fr(~ind,1),tone_fr(~ind,1),'ok','markerfacecolor','y')
            title('Baseline v. Tone Data Plot')
            xlabel('Baseline Firing Rates (spikes/second)')
            ylabel('Tone Task Firing Rate (spikes/second)')
            text(ul/20,.95*ul,['R = ' num2str(corr(base_fr(:,1),tone_fr(:,1)))])
        end
    
    % Tone v. ICMS data plot
    h2 = axes('units','pixels','position',pos(2,:));
        if isempty(tone_list) || isempty(icms_list)
            text(0.1,0.5,'Data not Available for Comparison')
            set(h2,'xtick',[],'xticklabel','','ytick',[],'yticklabel','')
            title('Tone v. ICMS Data Plot')
        else
            ul = ceil(max([tone_fr(:,1)+icms_fr(:,2); icms_fr(:,1)+tone_fr(:,2)]));
            set(h2,'xlim',[0 ul],'ylim',[0 ul])
            line([0 ul],[0 ul],'color','k','linewidth',2,'linestyle','--')
            line([tone_fr(:,1)'; tone_fr(:,1)'],[icms_fr(:,1)'-tone_fr(:,2)'; icms_fr(:,1)'+ tone_fr(:,2)'],'color','k','linewidth',2)
            line([tone_fr(:,1)'-icms_fr(:,2)'; tone_fr(:,1)'+ icms_fr(:,2)'],[icms_fr(:,1)'; icms_fr(:,1)'],'color','k','linewidth',2)
            hold on
            
            % Separate Units by channels 1-16 and 17-32
            ind = zeros(length(tone_list),1);
            for j = 1:length(tone_list)
                temp = str2double(tone_list{j}(5:6));
                if temp < 17
                    ind(j) = 1;
                end
            end
            ind = logical(ind);
            
            % Now plot it out
            plot(tone_fr(ind,1),icms_fr(ind,1),'ok','markerfacecolor','c'),
            plot(tone_fr(~ind,1),icms_fr(~ind,1),'ok','markerfacecolor','y')
            title('Tone v. ICMS Data Plot')
            xlabel('Tone Task Firing Rates (spikes/second)')
            ylabel('ICMS Task Firing Rate (spikes/second)')
            text(ul/20,.95*ul,['R = ' num2str(corr(tone_fr(:,1),icms_fr(:,1)))])
        end
    % Baseline v. ICMS data plot
    h3 = axes('units','pixels','position',pos(3,:));
        if isempty(base_list) || isempty(icms_list)
            text(0.1,0.5,'Data not Available for Comparison')
            set(h3,'xtick',[],'xticklabel','','ytick',[],'yticklabel','')
            title('Baseline v. ICMS Data Plot')
        else
            ul = ceil(max([base_fr(:,1)+icms_fr(:,2); icms_fr(:,1)+base_fr(:,2)]));
            set(h3,'xlim',[0 ul],'ylim',[0 ul])
            line([0 ul],[0 ul],'color','k','linewidth',2,'linestyle','--')
            line([base_fr(:,1)'; base_fr(:,1)'],[icms_fr(:,1)'-base_fr(:,2)'; icms_fr(:,1)'+ base_fr(:,2)'],'color','k','linewidth',2)
            line([base_fr(:,1)'-icms_fr(:,2)'; base_fr(:,1)'+ icms_fr(:,2)'],[icms_fr(:,1)'; icms_fr(:,1)'],'color','k','linewidth',2)
            hold on
            
            % Separate Units by channels 1-16 and 17-32
            ind = zeros(length(base_list),1);
            for j = 1:length(base_list)
                temp = str2double(base_list{j}(5:6));
                if temp < 17
                    ind(j) = 1;
                end
            end
            ind = logical(ind);
            
            % Now plot it out
            plot(base_fr(ind,1),icms_fr(ind,1),'ok','markerfacecolor','c'),
            plot(base_fr(~ind,1),icms_fr(~ind,1),'ok','markerfacecolor','y')
            title('Baseline v. ICMS Data Plot')
            xlabel('Baseline Firing Rates (spikes/second)')
            ylabel('ICMS Task Firing Rate (spikes/second)')
            text(ul/20,.95*ul,['R = ' num2str(corr(base_fr(:,1),icms_fr(:,1)))])
        end
        temp   = strfind(ICMS_files{q},'_');
        fileID = ICMS_files{q}(1:temp(2)-1);
        fileID(fileID == '_') = ' ';
        annotation('textbox',[.425 .85 .05 .05],'string',fileID,'fontsize',14,'edgecolor','none')
        drawnow
end
