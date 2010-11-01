function [pre_stim_psth,pre_trial_psth,unitIds] = change_pre_stim_fr(file,win)

% Requirements: this code requires the m-file 'data_matrix' to be in your
% path in order to run.
%
% Intent: To evaluate whether the pre-stimulus firing rates of the units changes
% prior to ICMS as the subjects learned the task. I'll compare the firing rates
% of the units at the trial start to the pre-stimulus period

% Inputs: 
% file  = a string providing the file of interest's name.
% bin   = bin in milliseconds
% win   = an integer number of ms prior to each event

if nargin < 1
    error('There must be at least one input: a string indicating a .mat file.')
end

if nargin < 3
    win  = 1000;
end

%%%%%%%%%%%%%%%%%
%%% Load file %%%
%%%%%%%%%%%%%%%%%
load(file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check for necessary inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Event005: Enter nosepoke event %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('Event015','var')
    error('The file does not contain the variable "Event015"')
end

% Register test event
event1 = Event015;

%%% Event011: Start trial event use for baseline %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('Event011','var')
    error('The file does not contain the variable "Event011"')
end

% A little preprocessing on Event011
% Account for double pulse Event011
ind = ones(1,length(Event011));
for n=2:length(Event011)
    if Event011(n) - Event011(n-1) < 0.8
        ind(n)   = 0;
    end
end

% Reigster baseline event
event2 = Event011(logical(ind));

% Grab Unit names
unitIds = who('sig*');

% Convert to timestamps to data matrix
% use data matrix to generate spikes
data_matrix
num_units = size(spikes,2);

%%% Generate PETH in response to nose poke %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find indices of event
epoch_ind = zeros(length(event1),1);
ind       = 0;

% Copy of time_line for speed up
time      = time_line;

for n = 1:length(event1)
    
    temp         = find(time>event1(n),1,'first');
    epoch_ind(n) = ind + temp;
    
    time         = time(temp+1:end);
    ind          = epoch_ind(n);
end

% Correct for terminal edge
if epoch_ind(end) + win/2 - 1 > size(spikes,1)
    % Toss trial out
    epoch_ind = epoch_ind(1:end-1);
end

% Generate peri-event response histogram for each unit    
% Grab Peri-Event Data
event_data = zeros(win,num_units,'uint16');
for q = 1:length(epoch_ind)
    event_data = event_data + uint16(spikes(epoch_ind(q)-win+1:epoch_ind(q),:));
end

% Convert average to spikes/second
pre_stim_psth   = 1000*single(event_data)./length(epoch_ind);

%%% Generate baseline unit firing rates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find indices of event
epoch_ind = zeros(length(event2),1);
ind       = 0;

for n = 1:length(event2)
    
    temp = find(time_line>event2(n),1,'first');
    epoch_ind(n) = ind + temp;
    
    time_line  = time_line(temp+1:end);
    ind        = epoch_ind(n);
end

% Generate pre-trial PETH histograms for each unit    
% Grab Peri-Event Data
event_data = zeros(win,num_units,'uint16');
for q = 1:length(epoch_ind)
    event_data = event_data + uint16(spikes(epoch_ind(q)-win/2:epoch_ind(q)+win/2-1,:));
end

% Convert average to spikes/second
pre_trial_psth = 1000*single(event_data)./length(epoch_ind);

