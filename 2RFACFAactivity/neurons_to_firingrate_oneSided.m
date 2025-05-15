function [FR,Nnew] = neurons_to_firingrate_oneSided(N)

% takes in the neurons data structure and outputs a firing rate matrix
% INPUT: N        - neurons data structre from the preprocessing output.
%        duration - the duration of the experiment, just put in the
%                   #columns of your EMG matrix or something

% OUTPUT: FR   - the firing rate matrix, still with zero rows
%         Nnew - the new neurons data structure with a new field ifr which
%                is the instantaneous firing rate vector (rows of FR)

addpath(genpath('Z:\Scripts\'));


units = 1:size(N,2);
trains = {N.train};
duration = max(cell2mat(trains'));

srate = 1000;                              % Hz 
min_timevec = 0;                           % sec
max_timevec = round(duration/1000);         % sec
sigma = 0.01;                              % sec (s.d. of Gaussian)
peak = 0;

Nnew = N;

columns = max_timevec*1000;
FR = zeros(length(N),columns);
pb = CmdLineProgressBar(' getting firing rate... ');
for i=1:length(N)
    pb.print(i, length(N));

    if(ismember(i,units))
        timestamps = trains{i}/1000; % convert from spike time in samples to seconds
        [fr,~,~] = smooth_spikes_oneSided(timestamps,srate,min_timevec,max_timevec,sigma,peak);
        fr = reshape(fr,1,[]);
        FR(i,:) = fr(1:columns);
    else
        fr = zeros(1,columns);
    end

    if(sum(fr) ~= 0)
        Nnew(i).ifr = fr(1:columns);
    end
end


disp('yay!');
% fprintf('--- \n Note: If you want to remove the rows that are all zero, just do FR( ~any(FR,2), : ) = []; \n Just remember that the unit indices become meaningless \n---');


end