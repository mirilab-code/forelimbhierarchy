function [bounds, climbing_bouts] = diya_get_inits(analogin, varargin)

%varargin pre

nInputs = numel(varargin);

if nInputs == 2
    pre = varargin{1};
    post = varargin{2};
elseif nInputs == 0
    pre = 0;
    post = 0;
else
    disp('something weird is inputted.')
    pre =0;
    post= 0;
end

set(0,'DefaultFigureVisible','off')
encoder_raw = analogin(2,:);


immobility_min_duration = 0.4; %minimal duration of "immobility" periods (sec)
immobility_min_duration_before_climbing = 1; %minimal duration of the immobility period preceding "climbing" bouts (sec)
min_climbing_angle = 10; %minimum angular distance travelled during a climbing bout (°)
[climbing_trials,encoder_filtered,isclimbing]= ...
        pre_processing_wheel(encoder_raw,immobility_min_duration,immobility_min_duration_before_climbing,min_climbing_angle);

init_times = ([climbing_trials.onset] + pre)';
if post==0
    stop_times = [climbing_trials.offset]';
else
    stop_times = ([climbing_trials.onset] + post)';
end

bounds = [init_times stop_times];

climbing_bouts = zeros(length(encoder_raw), 1, 'logical');
for i=1:length(init_times)
    climbing_bouts(init_times(i):stop_times(i)) = logical(1);
end































%%