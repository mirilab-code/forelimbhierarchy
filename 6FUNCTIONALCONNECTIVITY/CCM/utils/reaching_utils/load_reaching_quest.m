function [cfa_trains,rfa_trains, reaches] = load_reaching_quest(events_path, reaches_path)

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

