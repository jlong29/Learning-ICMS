%% SNR MODEL: WHISKER MODULATION: TEST HEMISPHERE
disp('SNR MODEL: WHISKER MODULATION: TEST HEMISPHERE')
% Process flow:
% 1) Designate subject
% 2) cd to subject's directory
% 3) Generate list of relevant files
% 4) Generate whisk_psth and mod_depth using 'whicker_mod_depth_PSTH.m'
% and after running all subjects save the data

%%%%%%%
% R29 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
file_names_R29 = {'R29_050808_spike&Events17_32.mat';
                  'R29_051008_spike&Events17_32.mat';
                  'R29_051208_spike&Events17_32.mat'};
              
for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_R29{n},20,1000,1000);
    eval(['R29_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['R29_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['R29_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['R29_mod_depth_i_' num2str(n) ' = mod_depth_i;'])

    clear whisk* mod* 
end
close all

%%%%%%
% V7 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
file_names_V7 = {'V7_032709_spikes&Events1_16.mat';
                 'V7_032909_spikes&Events1_16.mat';
                 'V7_033109_spikes&Events1_16.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V7{n},20,1000,1000);
    eval(['V7_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V7_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V7_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V7_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V8 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
file_names_V8 = {'V8_032909_spikes&Events1_16.mat';
                 'V8_033109_spikes&Events1_16.mat';
                 'V8_040109_spikes&Events1_16.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V8{n},20,1000,1000);
    eval(['V8_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V8_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V8_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V8_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V1 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
file_names_V1 = {'V1_092408_spikes&Events_1_16.mat';
                 'V1_092608_spikes&Events_1_16.mat';
                 'V1_093008_spikes&Events_1_16.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V1{n},20,1000,1000);
    eval(['V1_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V1_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V1_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V1_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V4 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
file_names_V1 = {'V4_092408_spikes&Events_1_16.mat';
                 'V4_092908_spikes&Events_1_16.mat';
                 'V4_100108_spikes&Events_1_16.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V1{n},20,1000,1000);
    eval(['V4_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V4_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V4_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V4_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
save('Figure5_whisk_mod_depth_test_bin20_win1000.mat','R29_*','V1_*','V4_*','V7_*','V8_*')
clear

%% SNR MODEL: WHISKER MODULATION: CONTROL HEMISPHERE
disp('SNR MODEL: WHISKER MODULATION: CONTROL HEMISPHERE')

% Process flow:
% 1) Designate subject
% 2) cd to subject's directory
% 3) Generate list of relevant files
% 4) Generate whisk_psth and mod_depth using 'whicker_mod_depth_PSTH.m'
% and after running all subjects save the data

%%%%%%%
% R29 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
file_names_R29 = {'R29_050808_spike&Events1_16.mat';
                  'R29_051008_spike&Events1_16.mat';
                  'R29_051208_spike&Events1_16.mat'};
              
for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_R29{n},20,1000,1000);
    eval(['R29_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['R29_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['R29_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['R29_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V7 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
file_names_V7 = {'V7_032709_spikes&Events17_32.mat';
                 'V7_032909_spikes&Events17_32.mat';
                 'V7_033109_spikes&Events17_32.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V7{n},20,1000,1000);
    eval(['V7_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V7_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V7_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V7_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V8 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
file_names_V8 = {'V8_032909_spikes&Events17_32.mat';
                 'V8_033109_spikes&Events17_32.mat';
                 'V8_040109_spikes&Events17_32.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V8{n},20,1000,1000);
    eval(['V8_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V8_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V8_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V8_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V1 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
file_names_V1 = {'V1_092408_spikes&Events_17_32.mat';
                 'V1_092608_spikes&Events_17_32.mat';
                 'V1_093008_spikes&Events_17_32.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V1{n},20,1000,1000);
    eval(['V1_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V1_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V1_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V1_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

%%%%%%
% V4 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
file_names_V1 = {'V4_092408_spikes&Events_17_32.mat';
                 'V4_092908_spikes&Events_17_32.mat';
                 'V4_100108_spikes&Events_17_32.mat'};

for n = 1:3
    [whisk_psth,mod_depth,mod_depth_e,mod_depth_i,ave_fr,unitIds] = ...
        whisker_mod_depth_PSTH(file_names_V1{n},20,1000,1000);
    eval(['V4_whisk_psth_' num2str(n) ' = whisk_psth;'])
    eval(['V4_mod_depth_t_' num2str(n) ' = mod_depth;'])
    eval(['V4_mod_depth_e_' num2str(n) ' = mod_depth_e;'])
    eval(['V4_mod_depth_i_' num2str(n) ' = mod_depth_i;'])
    
    clear whisk* mod*
end
close all

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig6_controls')
save('Figure6_whisk_mod_depth_control_bin20_win1000.mat','R29_*','V1_*','V4_*','V7_*','V8_*')
clear

%% SNR MODEL: INCREASED SIGNAL/DECREASED BACKGROUND
disp('SNR MODEL: INCREASED SIGNAL/DECREASED BACKGROUND')
% Process flow:
% 1) Designate subject
% 2) cd to subject's directory
% 3) Generate list of relevant files
% 4) Generate unit_psth to ICMS and using 'units_all_PSTH.m'
% 5) Generate SNR metrics
% and after running all subjects save the data

% Basic parameters
time_line = -500:10:490;
bin       = 10;
win       = 500;
ind       = length(time_line)/2;

%%%%%%%
% R29 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
file_names_R29 = {'R29_050808_spike&Events17_32.mat';
                  'R29_051008_spike&Events17_32.mat';
                  'R29_051208_spike&Events17_32.mat'};
              
for n = 1:3
    % Generate PSTHs for each unit reference to ICMS
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file_names_R29{n},'Event015',bin,win,win);
    eval(['R29_ICMS_psth_' num2str(n) ' = all_psth;'])
    
    % Calculate signal modulation
    data    = all_psth;
    
    % Determine size of variable and allocate memory
    temp    = size(data,2);
    
    %%%%%%%%%%%%%%%%%%%
    %%% SNR Metrics %%%
    %%%%%%%%%%%%%%%%%%%
    % Calculate SNR signal model
    sig_mod = zeros(temp,1);
    
    % Calculate SNR internal noise model
    nois_mod = zeros(temp,1);
    
    for p = 1:temp
        pre        = data(1:ind,p);
        post       = data(ind+1:end,p);
        
        % Calculate increase in signal relative to background (z-score)
        sig_mod(p) = (max(post) - mean(pre))/std(pre);
        
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
        nois_mod(p) = (mean(test)-mean(base))/...
                      (sqrt(var(test)/length(test) + var(base)/length(base)));
        
    end

    eval(['R29_ICMS_psth_' num2str(n) '_SNRs  = sig_mod;'])
    eval(['R29_ICMS_psth_' num2str(n) '_SNRb  = nois_mod;'])
    
end
close all

%%%%%%
% V7 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
file_names_V7 = {'V7_032709_spikes&Events1_16.mat';
                 'V7_032909_spikes&Events1_16.mat';
                 'V7_033109_spikes&Events1_16.mat'};

for n = 1:3
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file_names_V7{n},'Event015',bin,win,win);
    eval(['V7_ICMS_psth_' num2str(n) ' = all_psth;'])
    
    % Calculate signal modulation
    data    = all_psth;
    
    % Determine size of variable and allocate memory
    temp    = size(data,2);
    
    %%%%%%%%%%%%%%%%%%%
    %%% SNR Metrics %%%
    %%%%%%%%%%%%%%%%%%%
    % Calculate SNR signal model
    sig_mod = zeros(temp,1);
    
    % Calculate SNR internal noise model
    nois_mod = zeros(temp,1);
    
    for p = 1:temp
        pre        = data(1:ind,p);
        post       = data(ind+1:end,p);
        
        % Calculate increase in signal relative to background (z-score)
        sig_mod(p) = (max(post) - mean(pre))/std(pre);
        
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
        nois_mod(p) = (mean(test)-mean(base))/...
                      (sqrt(var(test)/length(test) + var(base)/length(base)));
        
    end

    eval(['V7_ICMS_psth_' num2str(n) '_SNRs  = sig_mod;'])
    eval(['V7_ICMS_psth_' num2str(n) '_SNRb  = nois_mod;'])
    
end
close all

%%%%%%
% V8 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
file_names_V8 = {'V8_032909_spikes&Events1_16.mat';
                 'V8_033109_spikes&Events1_16.mat';
                 'V8_040109_spikes&Events1_16.mat'};

for n = 1:3
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file_names_V8{n},'Event015',bin,win,win);
    eval(['V8_ICMS_psth_' num2str(n) ' = all_psth;'])
    
    % Calculate signal modulation
    data    = all_psth;
    
    % Determine size of variable and allocate memory
    temp    = size(data,2);
    
    %%%%%%%%%%%%%%%%%%%
    %%% SNR Metrics %%%
    %%%%%%%%%%%%%%%%%%%
    % Calculate SNR signal model
    sig_mod = zeros(temp,1);
    
    % Calculate SNR internal noise model
    nois_mod = zeros(temp,1);
    
    for p = 1:temp
        pre        = data(1:ind,p);
        post       = data(ind+1:end,p);
        
        % Calculate increase in signal relative to background (z-score)
        sig_mod(p) = (max(post) - mean(pre))/std(pre);
        
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
        nois_mod(p) = (mean(test)-mean(base))/...
                      (sqrt(var(test)/length(test) + var(base)/length(base)));
        
    end

    eval(['V8_ICMS_psth_' num2str(n) '_SNRs  = sig_mod;'])
    eval(['V8_ICMS_psth_' num2str(n) '_SNRb  = nois_mod;'])
    
end
close all

%%%%%%
% V1 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
file_names_V1 = {'V1_092408_spikes&Events_1_16.mat';
                 'V1_092608_spikes&Events_1_16.mat';
                 'V1_093008_spikes&Events_1_16.mat'};

for n = 1:3
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file_names_V1{n},'Event015',bin,win,win);
    eval(['V1_ICMS_psth_' num2str(n) ' = all_psth;'])
    
    % Calculate signal modulation
    data    = all_psth;
    
    % Determine size of variable and allocate memory
    temp    = size(data,2);
    
    %%%%%%%%%%%%%%%%%%%
    %%% SNR Metrics %%%
    %%%%%%%%%%%%%%%%%%%
    % Calculate SNR signal model
    sig_mod = zeros(temp,1);
    
    % Calculate SNR internal noise model
    nois_mod = zeros(temp,1);
    
    for p = 1:temp
        pre        = data(1:ind,p);
        post       = data(ind+1:end,p);
        
        % Calculate increase in signal relative to background (z-score)
        sig_mod(p) = (max(post) - mean(pre))/std(pre);
        
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
        nois_mod(p) = (mean(test)-mean(base))/...
                      (sqrt(var(test)/length(test) + var(base)/length(base)));
        
    end

    eval(['V1_ICMS_psth_' num2str(n) '_SNRs  = sig_mod;'])
    eval(['V1_ICMS_psth_' num2str(n) '_SNRb  = nois_mod;'])
    
end
close all

%%%%%%
% V4 %
%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
file_names_V1 = {'V4_092408_spikes&Events_1_16.mat';
                 'V4_092908_spikes&Events_1_16.mat';
                 'V4_100108_spikes&Events_1_16.mat'};

for n = 1:3
    [all_psth,ave_fr,unitIds,event_time] = units_all_PSTH(file_names_V1{n},'Event015',bin,win,win);
    eval(['V4_ICMS_psth_' num2str(n) ' = all_psth;'])
    
    % Calculate signal modulation
    data    = all_psth;
    
    % Determine size of variable and allocate memory
    temp    = size(data,2);
    
    %%%%%%%%%%%%%%%%%%%
    %%% SNR Metrics %%%
    %%%%%%%%%%%%%%%%%%%
    % Calculate SNR signal model
    sig_mod = zeros(temp,1);
    
    % Calculate SNR internal noise model
    nois_mod = zeros(temp,1);
    
    for p = 1:temp
        pre        = data(1:ind,p);
        post       = data(ind+1:end,p);
        
        % Calculate increase in signal relative to background (z-score)
        sig_mod(p) = (max(post) - mean(pre))/std(pre);
        
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
        nois_mod(p) = (mean(test)-mean(base))/...
                      (sqrt(var(test)/length(test) + var(base)/length(base)));
        
    end

    eval(['V4_ICMS_psth_' num2str(n) '_SNRs  = sig_mod;'])
    eval(['V4_ICMS_psth_' num2str(n) '_SNRb  = nois_mod;'])
    
end
close all

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
save('Figure5_ICMS_PSTH_SNR_bin10_win500.mat','R29_*','V1_*','V4_*','V7_*','V8_*')
clear

%% Process Group Averages
disp('Process Group Averages')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR Models: Increasing Signal and Decreasing Noise %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the relevant file generated in previous cell
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
file = 'Figure5_ICMS_PSTH_SNR_bin10_win500.mat';
load(file,'*_SNRs')

disp('Summary of SNR signal Increase Model')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of SNR signal Increase Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab per subject data
S1 = who('R29_ICMS*');
S2 = who('V1_ICMS*');
S3 = who('V7_ICMS*');
S4 = who('V8_ICMS*');
S5 = who('V4_ICMS*');

% Determine lower bound on units per session per animal
minind = zeros(3,5);

for p = 1:3
    eval(['minind(p,1) = ' 'length(' S1{p} ');'])
    eval(['minind(p,2) = ' 'length(' S2{p} ');'])
    eval(['minind(p,3) = ' 'length(' S3{p} ');'])
    eval(['minind(p,4) = ' 'length(' S4{p} ');'])
    eval(['minind(p,5) = ' 'length(' S5{p} ');'])
end

%%% Find ensemble size constraint for all subjects %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mins] = min(min(minind));

%%% Number of resample ensembles to be drawn %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N           = 500;
group_means = zeros(N,3);

for q = 1:N
    
    r1S1 = randsample(R29_ICMS_psth_1_SNRs,mins,true);
    r1S2 = randsample(V1_ICMS_psth_1_SNRs,mins,true);
    r1S3 = randsample(V7_ICMS_psth_1_SNRs,mins,true);
    r1S4 = randsample(V8_ICMS_psth_1_SNRs,mins,true);
    r1S5 = randsample(V4_ICMS_psth_1_SNRs,mins,true);
    
    r2S1 = randsample(R29_ICMS_psth_2_SNRs,mins,true);
    r2S2 = randsample(V1_ICMS_psth_2_SNRs,mins,true);
    r2S3 = randsample(V7_ICMS_psth_2_SNRs,mins,true);
    r2S4 = randsample(V8_ICMS_psth_2_SNRs,mins,true);
    r2S5 = randsample(V4_ICMS_psth_1_SNRs,mins,true);
    
    r3S1 = randsample(R29_ICMS_psth_3_SNRs,mins,true);
    r3S2 = randsample(V1_ICMS_psth_3_SNRs,mins,true);
    r3S3 = randsample(V7_ICMS_psth_3_SNRs,mins,true);
    r3S4 = randsample(V8_ICMS_psth_3_SNRs,mins,true);
    r3S5 = randsample(V4_ICMS_psth_1_SNRs,mins,true);
    
    group1 = [r1S1; r1S2; r1S3; r1S4; r1S5];
    group2 = [r2S1; r2S2; r2S3; r2S4; r2S5];
    group3 = [r3S1; r3S2; r3S3; r3S4; r3S5];
    
    group_means(q,1) = nanmean(group1);
    group_means(q,2) = nanmean(group2);
    group_means(q,3) = nanmean(group3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of Signal Increase Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndatSig = length(group1);
S       = group_means;
Smeans  = mean(S);
Sstds   = std(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear *_SNRs

disp('Summary of SNR Noise Decrease Model')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of SNR Noise Decrease Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(file,'*_SNRb')

% Number of resamples
N           = 500;
group_means = zeros(N,3);

for q = 1:N
    r1S1 = randsample(R29_ICMS_psth_1_SNRb,mins,true);
    r1S2 = randsample(V1_ICMS_psth_1_SNRb,mins,true);
    r1S3 = randsample(V7_ICMS_psth_1_SNRb,mins,true);
    r1S4 = randsample(V8_ICMS_psth_1_SNRb,mins,true);
    r1S5 = randsample(V4_ICMS_psth_1_SNRb,mins,true);
    
    r2S1 = randsample(R29_ICMS_psth_2_SNRb,mins,true);
    r2S2 = randsample(V1_ICMS_psth_2_SNRb,mins,true);
    r2S3 = randsample(V7_ICMS_psth_2_SNRb,mins,true);
    r2S4 = randsample(V8_ICMS_psth_2_SNRb,mins,true);
    r2S5 = randsample(V4_ICMS_psth_1_SNRb,mins,true);
    
    r3S1 = randsample(R29_ICMS_psth_3_SNRb,mins,true);
    r3S2 = randsample(V1_ICMS_psth_3_SNRb,mins,true);
    r3S3 = randsample(V7_ICMS_psth_3_SNRb,mins,true);
    r3S4 = randsample(V8_ICMS_psth_3_SNRb,mins,true);
    r3S5 = randsample(V4_ICMS_psth_1_SNRb,mins,true);
    
    group1 = [r1S1; r1S2; r1S3; r1S4; r1S5];
    group2 = [r2S1; r2S2; r2S3; r2S4; r2S5];
    group3 = [r3S1; r3S2; r3S3; r3S4; r3S5];

    group_means(q,1) = nanmean(group1);
    group_means(q,2) = nanmean(group2);
    group_means(q,3) = nanmean(group3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of internal noise model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndatNI   = length(group1);
NI       = group_means;
NImeans  = mean(NI);
NIstds   = std(NI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear *_SNRb

disp('Summary of External Noise Model')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of External Noise Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the relevant files
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
file = 'Figure5_whisk_mod_depth_test_bin20_win1000.mat';
load(file)
time_line = -1000:20:990;

% Grab per subject data
% total modulation
S1t = who('R29_mod_depth_t_*');
S2t = who('V1_mod_depth_t_*');
S3t = who('V7_mod_depth_t_*');
S4t = who('V8_mod_depth_t_*');
S5t = who('V4_mod_depth_t_*');

%excitation
S1e = who('R29_mod_depth_e_*');
S2e = who('V1_mod_depth_e_*');
S3e = who('V7_mod_depth_e_*');
S4e = who('V8_mod_depth_e_*');
S5e = who('V4_mod_depth_e_*');

% inhibition
S1i = who('R29_mod_depth_i_*');
S2i = who('V1_mod_depth_i_*');
S3i = who('V7_mod_depth_i_*');
S4i = who('V8_mod_depth_i_*');
S5i = who('V4_mod_depth_i_*');

% Bootstrap: sample with replacement units for each stage and subject given
% the minimum constraint calculated above.
% Details: For each stage and subject grab the minimum number of units with
% replacement, then calculate the per stage mean. Do this N times and use
% the SD of this as a measure of the SE for the group.

% Number of resamples
N             = 500;

% Let's look at absolute modulation as well as excitatory and inhibitory
% modulation
group_means   = zeros(N,3);
group_means_e = zeros(N,3);
group_means_i = zeros(N,3);

% Grab the variables
% total
S11t = eval(S1t{1}); S21t = eval(S2t{1}); S31t = eval(S3t{1}); S41t = eval(S4t{1}); S51t = eval(S5t{1});
S12t = eval(S1t{2}); S22t = eval(S2t{2}); S32t = eval(S3t{2}); S42t = eval(S4t{2}); S52t = eval(S5t{2});
S13t = eval(S1t{3}); S23t = eval(S2t{3}); S33t = eval(S3t{3}); S43t = eval(S4t{3}); S53t = eval(S5t{3});

% excitatory
S11e = eval(S1e{1}); S21e = eval(S2e{1}); S31e = eval(S3e{1}); S41e = eval(S4e{1}); S51e = eval(S5e{1});
S12e = eval(S1e{2}); S22e = eval(S2e{2}); S32e = eval(S3e{2}); S42e = eval(S4e{2}); S52e = eval(S5e{2});
S13e = eval(S1e{3}); S23e = eval(S2e{3}); S33e = eval(S3e{3}); S43e = eval(S4e{3}); S53e = eval(S5e{3});

% inhibitory
S11i = eval(S1i{1}); S21i = eval(S2i{1}); S31i = eval(S3i{1}); S41i = eval(S4i{1}); S51i = eval(S5i{1});
S12i = eval(S1i{2}); S22i = eval(S2i{2}); S32i = eval(S3i{2}); S42i = eval(S4i{2}); S52i = eval(S5i{2});
S13i = eval(S1i{3}); S23i = eval(S2i{3}); S33i = eval(S3i{3}); S43i = eval(S4i{3}); S53i = eval(S5i{3});

for q = 1:N
    % Sample absolute modulation and take average
    r1S1 = randsample(S11t,mins,true);
    r1S2 = randsample(S21t,mins,true);
    r1S3 = randsample(S31t,mins,true);
    r1S4 = randsample(S41t,mins,true);
    r1S5 = randsample(S51t,mins,true);
    
    r2S1 = randsample(S12t,mins,true);
    r2S2 = randsample(S22t,mins,true);
    r2S3 = randsample(S32t,mins,true);
    r2S4 = randsample(S42t,mins,true);
    r2S5 = randsample(S52t,mins,true);
    
    r3S1 = randsample(S13t,mins,true);
    r3S2 = randsample(S23t,mins,true);
    r3S3 = randsample(S33t,mins,true);
    r3S4 = randsample(S43t,mins,true);
    r3S5 = randsample(S53t,mins,true);
    
    group1 = [r1S1; r1S2; r1S3; r1S4; r1S5];
    group2 = [r2S1; r2S2; r2S3; r2S4; r2S5];
    group3 = [r3S1; r3S2; r3S3; r3S4; r3S5];

    group_means(q,1) = nanmean(group1);
    group_means(q,2) = nanmean(group2);
    group_means(q,3) = nanmean(group3);
    
    % sample excitatory modulation and take average
    r1S1 = randsample(S11e,mins,true);
    r1S2 = randsample(S21e,mins,true);
    r1S3 = randsample(S31e,mins,true);
    r1S4 = randsample(S41e,mins,true);
    r1S5 = randsample(S51e,mins,true);
    
    r2S1 = randsample(S12e,mins,true);
    r2S2 = randsample(S22e,mins,true);
    r2S3 = randsample(S32e,mins,true);
    r2S4 = randsample(S42e,mins,true);
    r2S5 = randsample(S52e,mins,true);
    
    r3S1 = randsample(S13e,mins,true);
    r3S2 = randsample(S23e,mins,true);
    r3S3 = randsample(S33e,mins,true);
    r3S4 = randsample(S43e,mins,true);
    r3S5 = randsample(S53e,mins,true);
    
    group1 = [r1S1; r1S2; r1S3; r1S4; r1S5];
    group2 = [r2S1; r2S2; r2S3; r2S4; r2S5];
    group3 = [r3S1; r3S2; r3S3; r3S4; r3S5];

    group_means_e(q,1) = nanmean(group1);
    group_means_e(q,2) = nanmean(group2);
    group_means_e(q,3) = nanmean(group3);
    
    % sample inhibitory modulation and take average
    r1S1 = randsample(S11i,mins,true);
    r1S2 = randsample(S21i,mins,true);
    r1S3 = randsample(S31i,mins,true);
    r1S4 = randsample(S41i,mins,true);
    r1S5 = randsample(S51i,mins,true);
    
    r2S1 = randsample(S12i,mins,true);
    r2S2 = randsample(S22i,mins,true);
    r2S3 = randsample(S32i,mins,true);
    r2S4 = randsample(S42i,mins,true);
    r2S5 = randsample(S52i,mins,true);
    
    r3S1 = randsample(S13i,mins,true);
    r3S2 = randsample(S23i,mins,true);
    r3S3 = randsample(S33i,mins,true);
    r3S4 = randsample(S43i,mins,true);
    r3S5 = randsample(S53i,mins,true);
    
    group1 = [r1S1; r1S2; r1S3; r1S4; r1S5];
    group2 = [r2S1; r2S2; r2S3; r2S4; r2S5];
    group3 = [r3S1; r3S2; r3S3; r3S4; r3S5];

    group_means_i(q,1) = nanmean(group1);
    group_means_i(q,2) = nanmean(group2);
    group_means_i(q,3) = nanmean(group3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summary of external noise model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Absolute modulation
ndatNE    = length(group1);
NE        = group_means;
NEmeans   = mean(NE);
NEstds    = std(NE);

% Excitatory modulation
NE_e      = group_means_e;
NEmeans_e = mean(NE_e);
NEstds_e  = std(NE_e);

% Inhibitory modulation
NE_i      = group_means_i;
NEmeans_i = mean(NE_i);
NEstds_i  = std(NE_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save relevant files
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
save('Learning_ICMS_SNR_models','S','Smeans','Sstds','NI','NImeans','NIstds',...
     'NE*','NEmeans*','NEstds*','ndatSig','ndatNI','ndatNE')

% system('%windir%\system32\rundll32.exe PowrProf.dll, SetSuspendState')

%%
% Load the relevant file generated in previous cell
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
file = 'Figure5_ICMS_PSTH_SNR_bin10_win500.mat';
load(file)

file = 'Figure5_whisk_mod_depth_test_bin20_win1000.mat';
load(file)

scrnsz = get(0,'screensize');
figure('Position',scrnsz);
subplot(5,3,1),plot(cumsum(sort(R29_ICMS_psth_1_SNRs)),'b'),hold on
    title('SNR: Increase Signal')
    ylabel('Sorted Cumsum')
subplot(5,3,4),plot(cumsum(sort(V1_ICMS_psth_1_SNRs)),'b'),hold on
subplot(5,3,7),plot(cumsum(sort(V7_ICMS_psth_1_SNRs)),'b'),hold on
subplot(5,3,10),plot(cumsum(sort(V8_ICMS_psth_1_SNRs)),'b'),hold on
subplot(5,3,13),plot(cumsum(sort(V4_ICMS_psth_1_SNRs)),'b'),hold on


subplot(5,3,1),plot(cumsum(sort(R29_ICMS_psth_2_SNRs)),'r')
subplot(5,3,4),plot(cumsum(sort(V1_ICMS_psth_2_SNRs)),'r')
subplot(5,3,7),plot(cumsum(sort(V7_ICMS_psth_2_SNRs)),'r')
subplot(5,3,10),plot(cumsum(sort(V8_ICMS_psth_2_SNRs)),'r')
subplot(5,3,13),plot(cumsum(sort(V4_ICMS_psth_2_SNRs)),'r')

subplot(5,3,1),plot(cumsum(sort(R29_ICMS_psth_3_SNRs)),'g')
subplot(5,3,4),plot(cumsum(sort(V1_ICMS_psth_3_SNRs)),'g')
subplot(5,3,7),plot(cumsum(sort(V7_ICMS_psth_3_SNRs)),'g')
subplot(5,3,10),plot(cumsum(sort(V8_ICMS_psth_3_SNRs)),'g')
subplot(5,3,13),plot(cumsum(sort(V4_ICMS_psth_3_SNRs)),'g')
    xlabel('Sorted Units')
    
subplot(5,3,2),plot(cumsum(sort(R29_ICMS_psth_1_SNRs)),'b'),hold on
    title('SNR: Decrease Noise')
    ylabel('Sorted Cumsum')
subplot(5,3,5),plot(cumsum(sort(V1_ICMS_psth_1_SNRb)),'b'),hold on
subplot(5,3,8),plot(cumsum(sort(V7_ICMS_psth_1_SNRb)),'b'),hold on
subplot(5,3,11),plot(cumsum(sort(V8_ICMS_psth_1_SNRb)),'b'),hold on
subplot(5,3,14),plot(cumsum(sort(V4_ICMS_psth_1_SNRb)),'b'),hold on

subplot(5,3,2),plot(cumsum(sort(R29_ICMS_psth_2_SNRb)),'r')
subplot(5,3,5),plot(cumsum(sort(V1_ICMS_psth_2_SNRb)),'r')
subplot(5,3,8),plot(cumsum(sort(V7_ICMS_psth_2_SNRb)),'r')
subplot(5,3,11),plot(cumsum(sort(V8_ICMS_psth_2_SNRb)),'r')
subplot(5,3,14),plot(cumsum(sort(V4_ICMS_psth_2_SNRb)),'r')

subplot(5,3,2),plot(cumsum(sort(R29_ICMS_psth_3_SNRb)),'g')
subplot(5,3,5),plot(cumsum(sort(V1_ICMS_psth_3_SNRb)),'g')
subplot(5,3,8),plot(cumsum(sort(V7_ICMS_psth_3_SNRb)),'g')
subplot(5,3,11),plot(cumsum(sort(V8_ICMS_psth_3_SNRb)),'g')
subplot(5,3,14),plot(cumsum(sort(V4_ICMS_psth_3_SNRb)),'g')
    xlabel('Sorted Units')
    
subplot(5,3,3),plot(cumsum(sort(R29_mod_depth_t_1)),'b'),hold on
    title('SNR: Whisker Response')
    ylabel('Sorted Cumsum')
subplot(5,3,6),plot(cumsum(sort(V1_mod_depth_t_1)),'b'),hold on
subplot(5,3,9),plot(cumsum(sort(V7_mod_depth_t_1)),'b'),hold on
subplot(5,3,12),plot(cumsum(sort(V8_mod_depth_t_1)),'b'),hold on
subplot(5,3,15),plot(cumsum(sort(V4_mod_depth_t_1)),'b'),hold on

subplot(5,3,3),plot(cumsum(sort(R29_mod_depth_t_2)),'r')
subplot(5,3,6),plot(cumsum(sort(V1_mod_depth_t_2)),'r')
subplot(5,3,9),plot(cumsum(sort(V7_mod_depth_t_2)),'r')
subplot(5,3,12),plot(cumsum(sort(V8_mod_depth_t_2)),'r')
subplot(5,3,15),plot(cumsum(sort(V4_mod_depth_t_2)),'r')

subplot(5,3,3),plot(cumsum(sort(R29_mod_depth_t_3)),'g')
subplot(5,3,6),plot(cumsum(sort(V1_mod_depth_t_3)),'g')
subplot(5,3,9),plot(cumsum(sort(V7_mod_depth_t_3)),'g')
subplot(5,3,12),plot(cumsum(sort(V8_mod_depth_t_3)),'g')
subplot(5,3,15),plot(cumsum(sort(V4_mod_depth_t_3)),'g')
    xlabel('Sorted Units')
    
set(findobj(gcf,'Type','line'),'linewidth',2)
axis(findobj(gcf,'Type','axes'),'tight','square')