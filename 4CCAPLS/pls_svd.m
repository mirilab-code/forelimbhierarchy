function [scoreX,scoreY,normedcovs,varcapX,varcapY] = pls_svd(X,Y)

% inputs: X and Y have to have the COLUMNS be neurons (or whatever
%         variables) and the ROWS to be time (or observations)



% center the data so that each neuron's activity has mean=0
X = X - mean(X,1);
Y = Y - mean(Y,1);

totalVar = var([X Y]);
summedVar = sum(totalVar);

crosscov = X' * Y;

[U,S,V] = svd(crosscov);
% S is the diagonal matrix of singular values which are covariances

% normalize the covariance by the sum of the variances for all the neurons in both regions
normedcovs = diag(S) / summedVar;






scoreX = X * U;
scoreY = Y * V;

varcapX = var(scoreX);
varcapY = var(scoreY);

varcapX = varcapX / sum(varcapX);
varcapY = varcapY / sum(varcapY);

end