function [Aconv,Bconv,Adiv,Bdiv] = convdiv(A,B)

% dim=1 looks at how many neuron is have signficant TE to neuron j: convergence
Aconv = mean(A,1);
Bconv = mean(B,1);

% dim=2 looks at how neuron i sends signficant TE to many neuron js: divergence
Adiv = mean(A,2);
Bdiv = mean(B,2);

end