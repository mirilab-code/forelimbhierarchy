


N = neurons{1};

% lets say we want to exlude all neurons above 2000um, all neurons with a
% firing rate below 1Hz and all narrow waveform neurons.
D = [N.depth];     % get all the depths
H = [N.firingrate];
W = [N.width];

% make a stack of logical arrays with your desired conditions
E = [D>2000; H<1; W<13];

exclude = sum(E,1)>0;

% now make a new neurons datastructure that removes all the units you don't want.
Nnew = N;
Nnew(exclude) = [];





