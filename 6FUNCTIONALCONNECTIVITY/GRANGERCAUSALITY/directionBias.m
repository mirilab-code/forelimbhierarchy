dates = ['09222020';'09232020';'09302020';'10012020';'10022020';...
    '10062020';'10242020';'10262020';'10272020';'10282020';...
    '11102020';'11122020'];
num_sessions = size(dates,1);
all_source_spikes = zeros(num_sessions,1);
all_target_spikes = zeros(num_sessions,1);
all_source_more = zeros(num_sessions,1);
all_target_more = zeros(num_sessions,1);
all_same = zeros(num_sessions,1);
for session = 1:num_sessions
    date = dates(session,:);
    str = sprintf('C:\\Users\\mirilab\\Documents\\Adam\\Granger_Quest\\%s',date);
    cd(str)
    str = sprintf('%s_Out.mat',date);
    load(str)
    psi = OutStruct.Psi2;
    num_neurons = size(SpikeTrains,1); 
    lin_ind = find(psi==1);
    [row,col] = ind2sub([num_neurons,num_neurons],lin_ind);
    summed_spikes = zeros(num_neurons,1);
    for ii = 1:num_neurons
        summed_spikes(ii) = length(find(SpikeTrains(ii,:,:)==1));
    end
    source_spikes = 0;
    target_spikes = 0;
    source_more = 0;
    target_more = 0;
    same = 0;
    for ii = 1:length(row)
        source_spikes = source_spikes + summed_spikes(col(ii));
        target_spikes = target_spikes + summed_spikes(row(ii));
        if summed_spikes(row(ii))>summed_spikes(col(ii))
            target_more = target_more+1;
        elseif summed_spikes(row(ii))<summed_spikes(col(ii))
            source_more = source_more+1;
        else
            same = same+1;
        end
    end
    all_source_spikes(session) = source_spikes;
    all_target_spikes(session) = target_spikes;
    all_target_more(session) = target_more;
    all_source_more(session) = source_more;
    all_same(session) = same;
end
cd('C:\Users\mirilab\Documents\Adam\Granger_Quest')
diff = all_source_spikes-all_target_spikes;
diff2 = all_source_more-all_target_more;
diff_pct = diff./(all_source_spikes+all_target_spikes);
figure;
scatter(1:num_sessions,diff_pct)
yline(0);
title('Aggregate Spikes: (+) -> more source spikes, (-) -> more target spikes')
xlabel('Sessions')
ylabel('Percent Difference')
figure;
scatter(1:num_sessions,diff2)
yline(0);
title('Per Connection Bias: (+) -> more source spikes, (-) -> more target spikes')
xlabel('Sessions')
ylabel('Number of Connections Difference')
%%
disp('source spikes')
disp(all_source_spikes')
disp('target spikes')
disp(all_target_spikes')
disp('greater source cnxs')
disp(all_source_more')
disp('greater target cnxs')
disp(all_target_more')

