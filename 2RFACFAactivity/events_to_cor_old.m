function [cortex,indices]= events_to_cor_old(num_probes,events, t_low, t_high)

%INPUTS events and t for threshold (for now, 1500)
%OUTPUTS CFA and RFA trains that can be run through JCCG
%need to check that u1 and u2 match with histogram
%t_low = 1500
%t_high = 0
cortex = cell(num_probes,1);
indices = cell(num_probes,1);

for ii = 1:num_probes
    e = events{ii};
    u = max(e(:,3))- t_low; 
    v = max(e(:,3))- t_high;

    cortical_neurons = find(e(:,3)>=u & e(:,3)<=v);
    st = e(cortical_neurons,2);
    ind = e(cortical_neurons,1);


    A = [ind'; st'];
    ecor=A';

    cortex{ii} = events_to_train(ecor);
    indices{ii} = unique(ind);
end

