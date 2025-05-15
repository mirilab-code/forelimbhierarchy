function [cfa_trains,rfa_trains, reaches] = load_reaching(preprocess_path)

date_very_temp = split(preprocess_path, '\');
date_temp = split(date_very_temp(end-1), '_');
date = date_temp{1,1};

events_path = sprintf('%s\\events.mat', preprocess_path);
reaches_path = sprintf('%s\\%s_1000_reach_bounds_edit.mat', preprocess_path, date);

events = load(events_path).events;
reaches = load(reaches_path).reach_bounds_edit;

CFA = cell2mat(events(1, 1));
RFA = cell2mat(events(2, 1));

CFA_train_temp = cort_events_to_train(CFA);
RFA_train_temp = cort_events_to_train(RFA);

cfa_trains = CFA_train_temp(~cellfun('isempty',CFA_train_temp));
rfa_trains = RFA_train_temp(~cellfun('isempty',RFA_train_temp));



%LOAD_REACHING Summary of this function goes here
%   Detailed explanation goes here

end

