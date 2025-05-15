%% plotting


figure;
scatter(scoreCFA(:,1),scoreRFA(:,1), '.k')

%%
figure;
corrplot = tiledlayout(2,2);
title(corrplot, 'CFA and RFA projects onto canonical variables');
xlabel(corrplot, 'time (ms)');
corrplot.TileSpacing = 'compact';

nexttile;
hold on
plot(score1(:,1));
plot(score2(:,1));
hold off
title('CV 1');
xticks(0:200:1600);

nexttile;
hold on
plot(score1(:,2));
plot(score2(:,2));
hold off
title('CV 2');
xticks(0:200:1600);

nexttile;
hold on
plot(score1(:,3));
plot(score2(:,3));
hold off
title('CV 3');
xticks(0:200:1600);

nexttile;
hold on
plot(score1(:,4));
plot(score2(:,4));
hold off
title('CV 4');
xticks(0:200:1600);


figure;
plot(r);
xlabel('r values');

%%
% split the pca score matrices into reaches and grasp and color code
T = 1:200;

figure;
hold on
for i=0:3
%     disp(i);
    window = T + (i*200);
    chunk_CFA = scoreCFA(window,1);
    reach_CFA = chunk_CFA(1:100);
    grasp_CFA = chunk_CFA(101:end);
    chunk_RFA = scoreRFA(window,1);
    reach_RFA = chunk_RFA(1:100);
    grasp_RFA = chunk_RFA(101:end);
    
    scatter(reach_CFA,reach_RFA, '.r');
    scatter(grasp_CFA,grasp_RFA, '.b');
end














disp('done~')
%%