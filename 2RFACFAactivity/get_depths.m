function depths = get_depths(path_to_kilosort)
addpath(genpath('Z:\Scripts\neuropixel\cortexlab-spikes'));
addpath('Z:\Scripts\neuropixel\npy-matlab');
cd(path_to_kilosort);

templates = readNPY('templates.npy');
whitening_mat = readNPY('whitening_mat.npy');
coords = readNPY('channel_positions.npy');
ycoords = coords(:,2);
spike_templates = readNPY('spike_templates.npy');
scaling_amps = readNPY('amplitudes.npy');

[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(templates, whitening_mat, ycoords, spike_templates, scaling_amps)

depths = spikeDepths;

end