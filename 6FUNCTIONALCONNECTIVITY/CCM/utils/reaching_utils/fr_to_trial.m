function output = fr_to_trial(fr,reaches)
%FR_TO_TRIAL Summary of this function goes here
%   Detailed explanation goes here
num_trials = size(reaches, 1);
%num_samples = 0;

output = cell(1, num_trials);
for i=1:num_trials
    %trial_size = (reaches(i,2) - reaches(i,1))+1;
    %num_samples = num_samples + trial_size;
    output{1, i} = fr(reaches(i,1):reaches(i,2));
end

