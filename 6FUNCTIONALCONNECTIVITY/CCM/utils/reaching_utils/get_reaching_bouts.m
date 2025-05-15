function [bounds, reaching_bouts] = get_reaching_bouts(reach_bounds, total_time, varargin)

%double check this diya

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

reach_bounds = sort(reach_bounds);

init_times = reach_bounds(:,1) + pre;
%if post==0
%stop_times = reach_bounds(:,2);
%else
stop_times = reach_bounds(:,1) + post; %i just modified this
%end

bounds = [init_times stop_times];

reaching_bouts = zeros(total_time, 1, 'logical');
for i=1:length(init_times)
    reaching_bouts(init_times(i):stop_times(i)) = logical(1);
end

end

