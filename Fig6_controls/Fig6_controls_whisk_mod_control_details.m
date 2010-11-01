%% Figure5 controls Details: The Control Hemisphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Whisker Modulation in CONTROL hemisphere %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig6_controls')
load('Figure6_whisk_mod_depth_control_bin20_win1000.mat')

% Determine lower bound on units per session per animal
minind = zeros(3,5);

% Grab per subject data
% total modulation
S1t = who('R29_mod_depth_t_*');
S2t = who('V1_mod_depth_t_*');
S3t = who('V7_mod_depth_t_*');
S4t = who('V8_mod_depth_t_*');
S5t = who('V4_mod_depth_t_*');

for p = 1:3
    eval(['minind(p,1) = ' 'length(' S1t{p} ');'])
    eval(['minind(p,2) = ' 'length(' S2t{p} ');'])
    eval(['minind(p,3) = ' 'length(' S3t{p} ');'])
    eval(['minind(p,4) = ' 'length(' S4t{p} ');'])
    eval(['minind(p,5) = ' 'length(' S5t{p} ');'])
end

[mins] = min(min(minind));

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
ndatNE      = length(group1);
NE_C        = group_means;
NE_Cmeans   = mean(NE_C);
NE_CstdsC   = std(NE_C);

% Excitatory modulation
NE_C_e      = group_means_e;
NE_Cmeans_e = mean(NE_C_e);
NEs_Ctds_e  = std(NE_C_e);

% Inhibitory modulation
NE_C_i      = group_means_i;
NEmeans_C_i = mean(NE_C_i);
NEstds_C_i  = std(NE_C_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
% append relevant variables to main file
save('Learning_ICMS_SNR_models.mat','NE_C*','-append')
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig6_controls')