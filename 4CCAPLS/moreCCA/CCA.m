%% analyze CCA

animal = 'ss2';
lag = 0;

load(sprintf('output/%s_%d.mat',animal,lag));

CFA = eval([animal '.CFA']);
RFA = eval([animal '.RFA']); 

CFA(isnan(CFA))=0;
RFA(isnan(RFA))=0;


% Remove zero rows
CFA(all(~CFA,2),:) = [];
RFA(all(~RFA,2),:) = [];

%% do PCA first

[coefCFA,scoreCFA,latentCFA] = pca(CFA');
[coefRFA,scoreRFA,latentRFA] = pca(RFA');

cvar_RFA = cumsum(latentCFA)./sum(latentCFA);
cvar_CFA = cumsum(latentRFA)./sum(latentRFA);

figure;
tiledlayout(2,1);
nexttile;
plot(cvar_CFA);
nexttile;
plot(cvar_RFA);

%%
pcadimCFA = find(cvar_CFA>0.95,1)
pcadimRFA = find(cvar_RFA>0.95,1)


[CV1,CV2,r,score1,score2] = canoncorr(scoreCFA(:,1:pcadimCFA),scoreRFA(:,1:pcadimRFA));

% % %%
% var1 = diag(cov(score1)) / trace(cov(scoreCFA(:,1:pcadimCFA)));
% var2 = diag(cov(score2)) / trace(cov(scoreRFA(:,1:pcadimRFA)));


%% orthonormal stuff and getting var capture

% but now, orthonormalize CV1 into orthCV1
orth1 = orth(CV1);
orth2 = orth(CV2);

Q1 = scoreCFA(:,1:pcadimCFA)*orth1;
Q2 = scoreRFA(:,1:pcadimRFA)*orth2;
varcap1 = flipud(diag(cov(Q1)) / trace(cov(scoreCFA)));
varcap2 = flipud(diag(cov(Q2)) / trace(cov(scoreRFA)));
sum(varcap1)
sum(varcap2)
hold on
plot(varcap1)
plot(varcap2)
hold off
%%
figure;
corrplot = tiledlayout(2,2);
title(corrplot, 'CFA and RFA projections onto canonical variables');
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
ylabel('r value');
xlabel('canonical dimension');

%% separate the reaches and grasps and spouts, this will be long

% IMPORTANT: run the script spout_EMG_pca.m to get the muscle PCAs
warning('run the script spout_EMG_pca.m to get the muscle PCAs if you havent!');

figure;
tiledlayout(5,8);

which_cvs = [1 2 5 10];
% which_cvs = [1 2 3 4];
% which_cvs = [5 6 7 8];
for cv=which_cvs
    disp(cv);
    ymin = min([score1(:,cv);score2(:,cv)]);
    ymax = max([score1(:,cv);score2(:,cv)]);
    for sec=1:8
        window = (200*(sec-1))+1:sec*200;      
        nexttile;
        hold on
        plot(-100:99,score1(window,cv));
        plot(-100:99,score2(window,cv));
        hold off
        if(cv==1)
            if(mod(sec,2))
                title(sprintf('reach spout %d', (sec+1)/2));
            else
                title(sprintf('grasp spout %d', (sec)/2));
            end
        end
        if(sec==1)
            s = sprintf('CV %d \n vcCFA %0.3f, vcRFA %0.3f, r=%0.3f',cv,varcap1(cv),varcap2(cv),r(cv));
            ylabel(s);
            
        end
        ylim([ymin ymax]);
    end 
end

for i=1:8
    nexttile;
    window = (200*(i-1))+1:i*200;    
    ymin = min(EMGscore(:,1));
    ymax = max(EMGscore(:,1));

    
    plot(-100:99,EMGscore(window,1));
    if(i==1)
        ylabel('PC1 of EMG');
    end
    ylim([ymin ymax]);
end


vc_CFARFA = [varcap1 varcap2];
save(sprintf('results/%s_varcap_CFARFA.mat',animal),'vc_CFARFA');



disp(animal);
disp('done!');