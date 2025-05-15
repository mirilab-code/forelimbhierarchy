function [rates] = train_binner(neurons, analogin, binwidth)
%TRAIN_BINNER Summary of this function goes here
%   Detailed explanation goes here
num_neurons = size(neurons, 2);
num_samples = size(analogin, 2);

rates = zeros(num_neurons, (num_samples / binwidth));

for i=1:num_neurons
    this_neuron = neurons(i);
    this_train = this_neuron.train;
    indices = ceil(this_train / binwidth);

    for k=1:size(indices)
        rates(i, indices(k)) = rates(i, indices(k)) + 1;
    end

end

rates = rates ./ (binwidth/1000);


