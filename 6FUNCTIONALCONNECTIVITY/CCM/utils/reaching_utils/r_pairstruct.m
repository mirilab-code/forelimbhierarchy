function pair_cell = r_pairstruct(CFA_pair, RFA_pair, reaches)
%REACHING_STRUCT Summary of this function goes here
%   Detailed explanation goes here

% output should be 1000 col cell of trials x time series

num_pairs = size(CFA_pair, 1);
pair_cell = cell(num_pairs, 2);
for i=1:num_pairs
    pair_cell{i, 1} = fr_to_trial(CFA_pair(i, :), reaches);
    pair_cell{i, 2} = fr_to_trial(RFA_pair(i,:), reaches);
end

