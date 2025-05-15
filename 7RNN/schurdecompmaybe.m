Q1 = symmetric.Ws(:,:,1);
Q2 = unconstrained.Ws(:,:,2);

d1 = eig(Q1);
d2 = eig(Q2);

figure;
hold on
plot(d1,'o')
plot(d2,'o')

%%
sym_decomp = eig_decomp(symmetric.Ws);
unc_decomp = eig_decomp(unconstrained.Ws);

hold on
plot(sym_decomp,'-k')
plot(unc_decomp,'-r')

function X = eig_decomp(stack)
    
    X = [];

    for i=1:size(stack,3)
        w = stack(:,:,i);
        d = abs(eig(w)');
        X = [X; d];
    end

    X = X';
end