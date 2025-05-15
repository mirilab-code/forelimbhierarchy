function trial_traces = get_trial_traces(emg,reach_bounds,lift_bounds,grasp_bounds)
    nReaches = size(reach_bounds,1);
    nMuscles = size(emg,1);
    trial_length = (lift_bounds(2)-lift_bounds(1))+(grasp_bounds(2)-grasp_bounds(1)) +2;
    trial_traces = zeros(nMuscles,trial_length,nReaches);
%     trial_traces = cell(nReaches,1);
    
    lift_1 = lift_bounds(1);
    lift_2 = lift_bounds(2);
    grasp_1 = grasp_bounds(1);
    grasp_2 = grasp_bounds(2);
    
    for i=1:nReaches
        lift = reach_bounds(i,1);
        grasp = reach_bounds(i,2);
        lift_trial = emg(:,lift+lift_1:lift+lift_2);
        grasp_trial = emg(:,grasp+grasp_1:grasp+grasp_2);
        trial_traces(:,:,i) = [lift_trial grasp_trial];
    end
    
end