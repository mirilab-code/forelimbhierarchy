%% do CCA


addpath(genpath('Z:code/'));
%%


matFiles = dir('output/*RFA_CONTROL.mat');  %do this for control
% matFiles = dir('output/*RFA.mat');  %do this for re al calculation
files = {matFiles.name};

for i = 1:length(files)
    fname = files{i};
    d = load(['output/' files{i}]);
    d.session = fname;
    
    nRFA = size(d.RFA_RG,1);
    nCFA = size(d.CFA_RG,1);

    fprintf('%s, %d in RFA, %d in CFA \n',fname,nRFA,nCFA);

    if(i==1)
        D = d;
    else
        D = [D d];
    end
end

D = aggregate_animals(D);

lags = D(1).lags;
num_PCs = 25;
lag0 = find(lags==0);

%%

weighted = [];
unweighted = [];
first_comp = [];
totalvarcap_RFACFAlag0 = [];
for i=1:length(D)
% for i=6:6
    disp(i);
    A = D(i);
    M = cca_all_lags(A.RFA_RG, A.CFA_RG, A.lags, num_PCs);

    weighted = [weighted; M.weighted];
    unweighted = [unweighted; M.unweighted];
    first_comp = [first_comp; M.first_comp];

    totalvarcap_RFACFAlag0 = [totalvarcap_RFACFAlag0; [sum(M.varcapRFA{lag0}) sum(M.varcapCFA{lag0})]];
end

%%



%%
figure;
hold on
plot(lags,weighted', '-k');




%%
mean_weighted = mean(weighted);
sem_weighted = std(weighted)/sqrt(size(weighted,1));

mean_unweighted = mean(unweighted);
sem_unweighted = std(unweighted)/sqrt(size(unweighted,1));

mean_comp1 = mean(first_comp);
sem_comp1 = std(first_comp)/sqrt(size(first_comp,1));





%%
figure;
subplot(2,1,1);
hold on
boundedline(lags,mean_weighted,sem_weighted, '-b','alpha')
boundedline(lags,mean_unweighted,sem_unweighted, '-r','alpha')
boundedline(lags,mean_comp1,sem_comp1, '-k','alpha')
% legend({'weighted','','unweighted','','first component'})
% xlim([-30 30]);
sgtitle(sprintf('number of PCs: %d',num_PCs));

subplot(2,1,2);
hold on
boundedline(lags,mean_weighted,sem_weighted, '-b','alpha')
boundedline(lags,mean_unweighted,sem_unweighted, '-r','alpha')
boundedline(lags,mean_comp1,sem_comp1, '-k','alpha')
legend({'weighted','','unweighted','','first component'})
% ylim([0.5 1])
xlim([-30 30]);



%%


























%%
disp('done!!');