%% Granger text parsing


dates = ['03272020';'03302020';'03312020';'04012020';'04022020';'04032020';...
    '09222020';'09232020';'09302020';'10012020';'10022020';'10052020';...
    '10062020';'10242020';'10262020';'10272020';'10282020';'11102020';...
    '11122020';'11152020';'11172020'];

num_dates = size(dates,1);

slurm_nums = ['6552927';'6605846';'6550312';'0000000';'0000000';'6605849';...
    '6424616';'6424628';'6459817';'6459818';'6459840';'0000000';'0000000';...
    '6459951';'6553016';'6459957';'0000000';'6956772';'6532663';'6525095';...
    '6328467'];
%pre slurms
% slurm_nums = ['1791551';'1765925';'1797437';'0000000';'0000000';'1797511';...
%     '1797878';'1798228';'1907439';'1907471';'1994616';'0000000';'0000000';...
%     '1994470';'1907509';'1907517';'0000000';'1907423';'1793938';'1794098';...
%     '1794520'];
no_slurms = [4,5,12,13,17];
slurms = 1:num_dates;
slurms(no_slurms) = [];

addpath(genpath('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains'));
neural_warnings_num = zeros(1,num_dates);
neural_warnings = cell(num_dates,1);
neural_pair_warnings = cell(num_dates,1);
for index = slurms
    date = dates(index,:);
    slurm_str = slurm_nums(index,:);
    str = ['slurm-',slurm_str,'.out'];
    fileID = fopen(str);
    slurmout = textscan(fileID,'%s');
    slurmout = slurmout{1};

    warnings = zeros(length(slurmout),1);
    for ii = 1:length(slurmout)
        if strcmp(slurmout{ii},'Warning:')
            warnings(ii) = 1;
        end
    end
    disp(length(find(warnings>0)))

    targets = zeros(length(slurmout),1);
    for ii = 1:length(slurmout)
        if strcmp(slurmout{ii},'target')
            targets(ii) = 1;
        end
    end
    disp(length(find(targets>0)))

    historys = zeros(length(slurmout),1);
    for ii = 1:length(slurmout)
        if strcmp(slurmout{ii},'history')
            historys(ii) = 1;
        end
    end
    disp(length(find(historys>0)))


    warning_inds = find(warnings>0);
    target_inds = find(targets>0);
    history_inds = find(historys>0);

    neuron_inds = zeros(length(warning_inds),1);
    neurons = cell(length(warning_inds),2);
    if ~isempty(warning_inds)
        for ii = 1:length(warning_inds)
            temp_ind = warning_inds(ii);
            target_check = (-1*target_inds)+temp_ind;
            target_check(target_check<0) = NaN;
            history_check = (-1*history_inds)+temp_ind;
            history_check(history_check<0) = NaN;
            a = min(target_check);
            b = min(history_check);
            if isempty(a)
                a = 10000000;
            end
            if isempty(b)
                b = 10000000;
            end
            if a<b
                [~,check_ind] = min(target_check);
                neuron_inds(ii) = target_inds(check_ind)+2;
                neurons{ii,1} = slurmout(neuron_inds(ii));
                neurons{ii,1}{1}(1) = [];
                neurons{ii,1} = str2num(neurons{ii,1}{1});
                neurons{ii,2} = slurmout(neuron_inds(ii)+5);
                neurons{ii,2}{1}(1) = [];
                neurons{ii,2} = str2num(neurons{ii,2}{1});
            else
                [~,check_ind] = min(history_check);
                neuron_inds(ii) = history_inds(check_ind)-2;
                neurons{ii,1} = slurmout(neuron_inds(ii));
                neurons{ii,1}{1}(1) = [];
                neurons{ii,1} = str2num(neurons{ii}{1});
                neurons{ii,2} = 0;
            end
        end
        neurons = cell2mat(neurons);
        hist_inds = find(neurons(:,2)==0);
        unique_neurons = unique(neurons(hist_inds,1));
        disp(length(unique_neurons))
        neurons(hist_inds,:) = [];
        unique_neuron_pairs = neurons;
        neural_warnings_num(index) = length(unique_neurons);
        neural_warnings{index} = unique_neurons;
        neural_pair_warnings{index} = unique_neuron_pairs;
    else
        disp('No Warnings!!')
    end
end

disp(neural_warnings_num)