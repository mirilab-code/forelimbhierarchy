function [reach,grasp] = one_spout_RG(M,reach_grasp,window)

reach_times = reach_grasp(:,1);
grasp_times = reach_grasp(:,2);

reach = get_trials_and_avg(M,reach_times,window);
grasp = get_trials_and_avg(M,grasp_times,window);



end