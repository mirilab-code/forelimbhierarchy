function [lzr,ctrl] = binmatrix_trials(neuron,laser,control,lzr_bw,ctrl_bw)

ntrials = length(laser);

lzr = zeros(ntrials,lzr_bw);
ctrl = zeros(ntrials,ctrl_bw);

for i=1:ntrials
    t = laser(i);
    w = zeros(1,lzr_bw);
    spikes = neuron(neuron>t & neuron<t+lzr_bw);
    spikes = ceil(spikes-t);
    w(spikes) = 1;
    lzr(i,:) = w;
    
    t = control(i);
    w = zeros(1,ctrl_bw);
    spikes = neuron(neuron>t & neuron<t+(ctrl_bw));
    spikes = ceil(spikes-t);
    w(spikes) = 1;
    ctrl(i,:) = w;
    
end


end