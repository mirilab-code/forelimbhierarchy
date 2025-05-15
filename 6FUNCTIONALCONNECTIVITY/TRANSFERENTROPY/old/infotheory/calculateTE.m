function [peakTE, CI, TEdelays] = calculateTE(stacked_bintrain,delay)

asdf = SparseToASDF(stacked_bintrain, 1);  % the last 1 means to use 1ms bin (doesn't matter in this case)


% tic
[peakTE, CI, TEdelays] = ASDFTE(asdf, 1:delay); % Now it has delay of 1ms to 30ms
% toc

% peakTE(i,j) is the peak transfer entropy from unit i to unit j

% removing self connections
peakTE = peakTE - diag(diag(peakTE));
CI = CI - diag(diag(CI));

end