function [events, spike_train, firing_rate] = readNPdata(path_to_ks_output)
% from the phy documentation:
%
% Initially, before running phy, the spike-cluster and spike-template 
% assignments are identical. If spike_clusters.npy does not exist, it is 
% automatically copied from spike_templates.npy. When modifying the 
% spike-cluster assignments in phy, only spike_clusters.npy is modified, 
% while spike_temp  lates.npy remains unchanged.

% note the the events and spike train are given in samples but the firing rate is in seconds

addpath('C:\Users\mirilab\Box\Miri Lab Shared\Scripts\npy-matlab');
addpath('C:\Users\mirilab\Box\Miri Lab Shared\Scripts');
addpath(path_to_ks_output);

spike_times = double(readNPY('spike_times.npy'));             % array of spike times, units are samples so they're taken at whatever hz the recording was. to get seconds divide by 30,000
spike_clusters = double(readNPY('spike_clusters.npy'));       % array of cluster indexes for each spike time

% the cluster ids start at 0... I'm wondering if I should just +1 all of them...
events_total = [spike_clusters spike_times];

KSlabels = tdfread('cluster_KSlabel.tsv','\t');     % cluster_KSLabel.tsv is from kilosort, if we fixed stuff with phy we need to look at cluster_info.tsv
good_ind = find(KSlabels.KSLabel(:,1)=='g');        % because there's a column for each letter, so good->g
good = KSlabels.cluster_id(good_ind);

% now remove the spikes from clusters that aren't good
f = find(ismember(events_total(:,1),good));
events = events_total(f,:);

% make spike train from events
u = unique(events(:,1));
spike_train = {};

for i=1:length(u)
    % get all the spike times for this unit
    unit_spike_times = find(events(:,1)==i);
    spike_train{i} = events(unit_spike_times,2);
end
spike_train = spike_train';
% spike_train = spike_train(~cellfun(@isempty, spike_train));       % this might not be what we want..?

disp('done making train!');

events_ms(:,2) = events(:,2)/30;          % the imec sampling rate, usually its 30kHz but check if it isnt
duration = max(events_ms(:,2));

% get activity over time
srate = 1000;                              % Hz 
min_timevec = 0;                           % sec
max_timevec = ceil(duration/1000);         % sec
sigma = 0.01;                              % sec
peak = 0;

FR = zeros(length(u),max_timevec*1000+1);

for i=1:length(spike_train)
    disp(i);
    spkvec = spike_train{i}/1000;
    [fr,~,~] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
    FR(i,:) = fr';
end

firing_rate = FR;

disp('done getting activity!');


end