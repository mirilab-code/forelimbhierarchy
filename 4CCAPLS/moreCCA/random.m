[a,b,c] = pca(CFA');
v = c/sum(c);

%%
hold on
plot(cumsum(v))
yline(0.9)


%%
%%
excludeCFA = ~ismember(events{1}(:,1),CFA_units);
events{1}(excludeCFA,:) = [];



















%%