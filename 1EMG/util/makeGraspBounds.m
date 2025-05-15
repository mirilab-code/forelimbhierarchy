function grasp_bounds = makeGraspBounds(num_spouts,grasp_window,reach_bounds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
grasp_bounds = cell(1,num_spouts);
for ii = 1:num_spouts
   grasp_bounds{ii} = [reach_bounds{ii}(:,2)+grasp_window(1),...
       reach_bounds{ii}(:,2)+grasp_window(2)];
end

