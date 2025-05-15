function [shift_CFA_trial,shift_RFA_trial] = circshift_pair(CFA_FR, RFA_FR, reaching_bouts)

CFA_climbing= CFA_FR(reaching_bouts);
RFA_climbing = RFA_FR(reaching_bouts);

%this is incredibly stupid

temp = [diff(reaching_bouts);0]; % this should give a -1 when stitching
stitch_times_temp = find(temp(reaching_bouts)==-1);
stitch_times = [0; stitch_times_temp];

num_ts = length(CFA_climbing);
shifts = 3000:num_ts-3000;   % acceptable shifts are more than 3seconds or less than T-3 seconds
r = randi([1 2]); % randomly shift only one region
unit_shift = randsample(shifts,1);    % true = with replacement

if r == 1
    CFA_climbing = circshift(CFA_climbing, unit_shift);
    %fprintf('CFA shifted %d amount \n', unit_shift)
else
    RFA_climbing = circshift(RFA_climbing, unit_shift);
    %fprintf('RFA shifted %d amount \n', unit_shift)
end

num_trials = length(stitch_times) - 1;

shift_CFA_trial = cell(1, num_trials);
shift_RFA_trial = cell(1, num_trials);

for i=1:num_trials
    shift_CFA_trial{1,i} = CFA_climbing(stitch_times(i)+1:stitch_times(i+1));
    shift_RFA_trial{1,i} = RFA_climbing(stitch_times(i)+1:stitch_times(i+1));
end

