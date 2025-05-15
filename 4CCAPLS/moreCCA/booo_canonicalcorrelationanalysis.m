%% canonical correlation analysis

datapath = 'Z:\akiko\10242020_MA1_g0\preprocess_with_acg';
addpath('C:\Users\mirilab\Documents\Mark\util\');
addpath(datapath);

files = fullfile(datapath, '*.mat');
matFiles = dir(files);
for i = 1:length(matFiles)
    baseFileName = fullfile(datapath, matFiles(i).name);
    load(baseFileName);
end
disp('done loading!');
nProbes = length(events);



trains = {};
FR = {};
for i=1:nProbes
    trains{i} = events_to_train(events{i});
    FR{i} = train_to_FR(trains{i});
end
nCells1 = length(trains{1});
nCells2 = length(trains{2});
duration = length(FR{1}(1,:));
disp('done!');

%%
[CV1,CV2,r,score1,score2] = canoncorr(FR{1}',FR{2}');
%%

figure;
t = tiledlayout(2,2);
title(t,'Canonical Scores of X vs Canonical Scores of Y')
xlabel(t,'Canonical Variables of X')
ylabel(t,'Canonical Variables of Y')
t.TileSpacing = 'compact';

T = 1:1000;

nexttile
plot(score1(:,1),score2(:,1),'.')
xlabel('u1')
ylabel('v1')
title(sprintf('r = %f',r(1)));

nexttile
plot(score1(:,2),score2(:,2),'.')
xlabel('u2')
ylabel('v2')
title(sprintf('r = %f',r(2)));

nexttile
plot(score1(:,3),score2(:,3),'.')
xlabel('u3')
ylabel('v3')
title(sprintf('r = %f',r(3)));

nexttile
plot(score1(:,4),score2(:,4),'.')
xlabel('u4')
ylabel('v4')
title(sprintf('r = %f',r(4)));


%% shuffle FR2 and do CCA to see how significant the CVs are
shiftFR2 = [];
for i=1:size(FR{2},1)
    shift = randsample(1:duration,1);
    shiftFR2 = [shiftFR2; circshift(FR{2}(i,:),shift)];
end
%%

[CV1_shifted,CV2_shifted,r_shifted,score1_shifted,score2_shifted] = canoncorr(FR{1}',shiftFR2');

%%
figure;
hold on
histogram(r)
histogram(r_shifted);
hold off








disp('done!');
%%