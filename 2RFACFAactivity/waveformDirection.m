function waveform_dir = waveformDirection(channelMapFile,imec_path,ks_path,units)
%waveformDireaction Determines whether the extreme value of the template
%wavefrom for a given series of units is positive (0) or negative (1).
%   Detailed explanation goes here
ld = load(channelMapFile);
chanMap = Neuropixel.ChannelMap(channelMapFile);
imec = Neuropixel.ImecDataset(imec_path,'channelMap', chanMap);
ks = Neuropixel.KiloSortDataset(ks_path,'imecDataset',imec');
ks.load();
metrics = ks.computeMetrics();
unit_waveforms = metrics.cluster_waveform(units,:);
waveform_dir = zeros(1,length(units));
for ii = 1:size(unit_waveforms,1)
    waveform_sign = sign(unit_waveforms(ii,:));
    peak_index = find(abs(unit_waveforms(ii,:))==...
        max(abs(unit_waveforms(ii,:))));
    if waveform_sign(peak_index) == -1
        waveform_dir(ii) = 1;
    else
        waveform_dir(ii) = 0;
    end
end

