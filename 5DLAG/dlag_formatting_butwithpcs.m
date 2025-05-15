
clear all;

%% Load files
session_path = uigetdir('', 'Data Path');
cd(session_path)
date = session_path(end-7:end);

data = uigetfile('','PCs');
load(data);



%% Struct creation

numTrials = size(onset_bounds,1);

pc_struct = struct;
for t = 1:numTrials
    pc_struct(t).trialId = t;
    pc_struct(t).T = pre+post+1;
    pc_struct(t).spikes = [reshape(CFA_pcs(t,:,:), [size(CFA_pcs,3), size(CFA_pcs,2)]); reshape(RFA_pcs(t,:,:), [size(RFA_pcs,3), size(RFA_pcs,2)])];
end



%% getSeq

cd('Z:\Sarah\DLAG\')
addpath('DLAG-1.0.0\')

all_seq = getSeq(pc_struct, 1);


%% Save
save_path = uigetdir('', 'Save Location');
cd(save_path)

save(sprintf('%s_binned_spike_counts.mat', date), 'all_seq', 'CFA_excluded', 'RFA_excluded');


%% Notes

% Need: xDims_across and xDims_within from FA?

