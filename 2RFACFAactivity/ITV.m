function S = ITV(M)
% measure the InterTrial Variability. M is a 3d matrix (muscles x time x
% trials). IVT takes the standard deviation of a muscle at each time point
% and then takes the average of those stdevs for each muscle

num_muscles = size(M,1);
duration = size(M,2);
num_trials = size(M,3);

stdevs = zeros(num_muscles,duration);

for m=1:num_muscles
    for t=1:duration
        sigma = std(M(m,t,:));
        stdevs(m,t) = sigma;
    end
end 

% disp(size(stdevs));
S = mean(stdevs,2);


end