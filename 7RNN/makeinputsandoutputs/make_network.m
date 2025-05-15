%% make network architecture

% neuron types (all mutuallly exclusive):
% - Intercortical pyramidal (IC): 3% of each region, excitatory connections from one region to the other
% - Striatal projection pyramidal (SPN): 10% of each region, these are the output that makes the EMG trace
% - Input pyramidal (IN): 10% of each region, these are what we feed the go cue into


N_RFA = 500;
N_CFA = 500;
N_total = N_RFA + N_CFA;
perc_exc_RFA = 0.8;
perc_exc_CFA = 0.8;
N_exc_RFA = N_RFA * perc_exc_RFA;
N_exc_CFA = N_CFA * perc_exc_CFA;
N_inh_RFA = N_RFA - N_exc_RFA;
N_inh_CFA = N_CFA - N_exc_CFA;

broad_types = containers.Map([1 -1],{'exc','inh'});
flags_exc_RFA = zeros(N_exc_RFA,1)+1;
flags_inh_RFA = zeros(N_inh_RFA,1)-1;
flags.broad_RFA = [flags_exc_RFA; flags_inh_RFA];

flags_exc_CFA = zeros(N_exc_CFA,1)+1;
flags_inh_CFA = zeros(N_inh_CFA,1)-1;
flags.broad_CFA = [flags_exc_CFA; flags_inh_CFA];

%%
types = containers.Map([1 2 3 0],{'input','output','inter','not important'});

perc_input_RFA = 0.1;
perc_input_CFA = 0.1;
perc_output_RFA = 0.1;
perc_output_CFA = 0.1;

N_input_RFA = N_RFA * perc_exc_RFA * perc_input_RFA;
N_input_CFA = N_CFA * perc_exc_CFA * perc_input_CFA;
N_output_RFA = N_RFA * perc_exc_RFA * perc_output_RFA;
N_output_CFA = N_CFA * perc_exc_CFA * perc_output_CFA;
% N_inter_RFA = N_RFA * perc_exc_RFA * perc_inter_RFA;
% N_inter_CFA = N_CFA * perc_exc_CFA * perc_inter_CFA;

flags_input_RFA = zeros(N_input_RFA,1)+1;
flags_output_RFA = zeros(N_output_RFA,1)+2;
% flags_inter_RFA = zeros(N_inter_RFA,1)+3;
flags_therest_RFA = zeros(N_RFA - (N_input_RFA+N_output_RFA),1);

flags_input_CFA = zeros(N_input_CFA,1)+1;
flags_output_CFA = zeros(N_output_CFA,1)+2;
% flags_inter_CFA = zeros(N_inter_CFA,1)+3;
flags_therest_CFA = zeros(N_CFA - (N_input_CFA+N_output_CFA),1);


flags.RFA = [flags_input_RFA; flags_output_RFA; flags_therest_RFA];
flags.CFA = [flags_input_CFA; flags_output_CFA; flags_therest_CFA];

%%
popCFA = struct;
for i=1:N_CFA
    popCFA(i).region = 'cfa';
    v = flags.broad_CFA(i);
    btype = broad_types(v);
    popCFA(i).ei = btype;

    v = flags.CFA(i);
    t = types(v);
    popCFA(i).type = t;
end

popRFA = struct;
for i=1:N_CFA
    popRFA(i).region = 'rfa';
    v = flags.broad_RFA(i);
    btype = broad_types(v);
    popRFA(i).ei = btype;

    v = flags.RFA(i);
    t = types(v);
    popRFA(i).type = t;
end

pop = [popCFA popRFA]';

%%
% inputs: the population struct, 'cfa' or 'rfa', 'input' or 'output'
get_indices = @(P,region,type) ...
    strcmp({P.region},region) & strcmp({P.type},type);
get_ei = @(P,region,ei) ...
    strcmp({P.region},region) & strcmp({P.ei},ei);


D.pop = pop;
D.get_indices = get_indices;
D.get_ei = get_ei;
% now you can query the indices of certain neurons like cfainputs = D.get_indices(D.pop,'cfa','input')

%% make weight matrix
perc_inter_RFA = 0.03;
perc_inter_CFA = 0.03;

inputscfa = find(D.get_indices(D.pop,'cfa','input'));
outputscfa = find(D.get_indices(D.pop,'cfa','output'));
inputsrfa = find(D.get_indices(D.pop,'rfa','input'));
outputsrfa = find(D.get_indices(D.pop,'rfa','output'));
exc_cfa = find(D.get_ei(D.pop,'cfa','exc'));
exc_rfa = find(D.get_ei(D.pop,'rfa','exc'));
inh_cfa = find(D.get_ei(D.pop,'cfa','inh'));
inh_rfa = find(D.get_ei(D.pop,'rfa','inh'));
possible_inter_cfa = intersect(find(D.get_indices(D.pop,'cfa','not important')), exc_cfa);
possible_inter_rfa = intersect(find(D.get_indices(D.pop,'rfa','not important')), exc_rfa);

exc = union(exc_cfa,exc_rfa);
inh = union(inh_cfa,inh_rfa);
possible_inter = union(possible_inter_cfa, possible_inter_rfa)';

%%
intra_scale = 2;
inter_scale = 2;
W = zeros(N_total,N_total);

p_inter = perc_inter_CFA * perc_inter_RFA;
p_intra = 0.05;
% p_inter = 1;
for i=1:N_total
    for j=1:N_total
        i_region = D.pop(i).region;
        j_region = D.pop(j).region;
        
        % first check to do the intraregional synapse
        if(strcmp(i_region,j_region))
%             fprintf('%s, %s \n',i_region,j_region);
            i_ei = D.pop(i).ei;
            x = rand();
            if(strcmp(i_ei,'exc') && x <= p_intra)
                W(i,j) = abs(lognrnd(0,1));
            elseif(strcmp(i_ei,'inh') && x <= p_intra)
                W(i,j) = -intra_scale*abs(lognrnd(0,1));
            end
        end
        
        % then check to do interregional, but make sure no input or output neurons are connected.
        if(~strcmp(i_region,j_region))
            if(ismember(i,possible_inter))  % the source neuron has to be exc (just not input or output) but can synapse to exc or inh. If you want to make it so that these also only synapse onto non-input or output neurons, just add &&ismember(i,possible_inter) in the if condition.
                x = rand();
                if(x <= p_inter)
                    W(i,j) = inter_scale*abs(lognrnd(0,1));
                end
            end
        end

    end
end

D.W = W;

imagesc(W)

%%
network = D;
save('model_setup_data\network.mat','network');
disp('done saving!')












