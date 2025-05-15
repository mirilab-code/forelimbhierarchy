%% load data
load('trains.mat');
%%
which_train = 2;
st1 = trains{which_train};
trainlen = length(trains{which_train});
duration = max(max(cell2mat(trains{1}),max(cell2mat(trains{2}))));

% calculate firing rate of each unit
nspikes = cellfun(@length,st1);
firingrates = nspikes./duration;

% we want to get rid of these low firing units, not enough data for ISI/info
lowFR = nspikes < 2000;

st = st1(~lowFR);
% now we recalculate the firing rates 
nspikes = cellfun(@length,st);
firingrates = nspikes./duration;
N = length(st);

bin_width = 1;

% get ISIs and stuff
ISI = cell(N,1);
info = cell(N,1);
probdists = cell(N,1);
distparams = zeros(N,2);

for i=1:N
%     disp([i N]);
    train = st{i};
    if(train(end) ~= duration)
        train = [train; duration];
    end
    
    isi = diff(train);
    maxISI = max(isi);
    % get rid of any 0 isi's because that messes up logarithms
    isi0s = find(isi == 0);
    minnotzero = min(isi(isi>0));
    epsilonisi = minnotzero/10;   % just make any 0 isi's really small
    isi(isi0s) = epsilonisi;
    
    ISI{i} = isi;
    pd = fitdist(isi,'Gamma');
    probdists{i} = pd;
    a = pd.a;
    b = pd.b;
    params = [a b];
    distparams(i,:) = params;
    x = 0:maxISI;
    y = pdf(pd,x);
    y(y==0) = NaN;
    y(y==Inf) = NaN;
    yfill = fillmissing(y,'nearest');
    info{i} = -log2(yfill);
%     edges = 1:bin_width:max(isi)+1;
%     v = histcounts(isi,edges,'Normalization','probability');
%     sum(v);
%     % this prevents the info from going infinite
%     v(v==0) = NaN;
%     v = fillmissing(v,'nearest');
%     probdists{i} = v;
%     info{i} = -log2(v);
end

disp('done!');
%% get running information
running_info = cell(N,1);
tsls = cell(N,1);
for i=1:N
    disp(i);
    running_info{i} = info_train(ISI{i},info{i},probdists{i});
    tsls{i} = time_since_last_spike(ISI{i});
end
disp('done getting running info!');

%% get the information of each spike
info_per_spike = cell(N,1);
for i=1:N
    disp(i);
    info_per_spike{i} = spike_info(ISI{i},probdists{i});
end
%
% now to get bits per second we do bits/spike * spikes/second
bits_per_spike = cellfun(@sum,info_per_spike) ./ cellfun(@length,info_per_spike);
bits_per_second = bits_per_spike .* firingrates;
disp('done!');





















%%


















%%