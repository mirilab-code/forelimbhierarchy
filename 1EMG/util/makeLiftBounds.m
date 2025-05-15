function lift_bounds = makeLiftBounds(num_spouts,lift_window,reach_bounds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lift_bounds = cell(1,num_spouts);
for ii = 1:num_spouts
   lift_bounds{ii} = [reach_bounds{ii}(:,1)+lift_window(1),...
       reach_bounds{ii}(:,1)+lift_window(2)];
end

