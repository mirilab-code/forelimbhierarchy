function events = get_events_from_ks(path_to_ks,whichunits)
% INPUT: path_to_ks: the path to the kilosort output. It should be a folder with 
%        stuff like amplitudes.npy, channel_map.npy, etc
%        whichunits: 'good' or 'good and mua'. Kilosort labels each
%        cluster as 'good' or 'mua' so pick which ones you want to include.
% OUTPUT: [a lot]x3 matrix with the columns being cluster id, spike time, depth
% 
% from the phy documentation:
% Initially, before running phy, the spike-cluster and spike-template 
% assignments are identical. If spike_clusters.npy does not exist, it is 
% automatically copied from spike_templates.npy. When modifying the 
% spike-cluster assignments in phy, only spike_clusters.npy is modified, 
% while spike_templates.npy remains unchanged.

% note the the events are given in samples.
% With kilosort the units are microns and the 
%  channels with depth=0 are the ones nearest the tip.


% outputs: 
% - spikeAmps is length nSpikes vector with amplitude in unwhitened space
% of every spike
% - spikeDepths is the position along the probe of every spike (according
% to the position of the template it was extracted with)
% - templateDepths is the position along the probe of every template
% - templateAmps is the amplitude of each template
% - tempsUnW are the unwhitened templates
% - templateDuration is the trough-to-peak time (in samples)
% - waveforms: returns the waveform from the max-amplitude channel
%
% inputs: 
% - temps, the templates (nTemplates x nTimePoints x nChannels)
% - winv, the whitening matrix (nCh x nCh)
% - ycoords, the coordinates of the channels (nCh x 1)
% - spikeTemplates, which template each spike came from (nSpikes x 1)
% - tempScalingAmps, the amount by which the template was scaled to extract
% each spike (nSpikes x 1)

addpath('Z:\Scripts\neuropixel\npy-matlab');
addpath(path_to_ks);

spike_times = double(readNPY('spike_times.npy'));             % array of spike times, units are samples so they're taken at whatever hz the recording was. to get seconds divide by 30,000
spike_clusters = double(readNPY('spike_clusters.npy'));       % array of cluster indexes for each spike time

templates = readNPY('templates.npy');
whitening_mat = readNPY('whitening_mat.npy');
coords = readNPY('channel_positions.npy');
ycoords = coords(:,2);
spike_templates = readNPY('spike_templates.npy');
scaling_amps = readNPY('amplitudes.npy');

[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(templates, whitening_mat, ycoords, spike_templates, scaling_amps);

% disp([length(spike_clusters) length(spike_times) length(spikeDepths)]);


events_total = [spike_clusters spike_times spikeDepths];

KSlabels = tdfread('cluster_KSlabel.tsv','\t');     % cluster_KSLabel.tsv is from kilosort, if we fixed stuff with phy we need to look at cluster_info.tsv
if (isequal(whichunits, 'good'))
    good_ind = find(KSlabels.KSLabel(:,1)=='g');        % because there's a column for each letter, so good->g
    good = KSlabels.cluster_id(good_ind);
    % now remove the spikes from clusters that aren't good
    f = find(ismember(events_total(:,1),good));
    events = events_total(f,:);

elseif (isequal(whichunits, 'good and mua'))
    good_ind = find(KSlabels.KSLabel(:,1)=='g');        % because there's a column for each letter, so good->g
    good = KSlabels.cluster_id(good_ind);
    mua_ind = find(KSlabels.KSLabel(:,1)=='m');
    mua = KSlabels.cluster_id(mua_ind);
    inds = union(good,mua);
    % now remove the spikes from clusters that aren't good or mua
    f = find(ismember(events_total(:,1),inds));
    events = events_total(f,:);
    
elseif (isequal(whichunits, 'all'))
    events = events_total;
else
    error('second argument must be good, good and mua, or all')
end



% make the unit index start at 1 instead of 0
events(:,1) = events(:,1)+1;

end