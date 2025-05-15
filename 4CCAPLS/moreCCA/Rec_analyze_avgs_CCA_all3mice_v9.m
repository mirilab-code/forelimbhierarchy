%Look at Neural and muscle PCs or SVDs
%check that muscles have the same phase
close all
clear all
numFiles = 11;

fr_th = 1;
PCAdim = 4;
plotDim = 4;
smoWin = 100;
smoStd = 10; %should yield 10
numIters = 1000;

animals = {'C20','C31','C32'};
suffix = {'_10ms_500_abby','_10ms_500_abby','_10ms_500_abby'};

CorrsLoco = zeros(PCAdim,length(animals));
CorrsReach = zeros(PCAdim,length(animals));
NeuralVarCap_Reach_CVs = zeros(PCAdim,length(animals));
MuscVarCap_Reach_CVs = zeros(PCAdim,length(animals));
NeuralVarCap_Loco_CVs = zeros(PCAdim,length(animals));
MuscVarCap_Loco_CVs = zeros(PCAdim,length(animals));

for n = 1:3

    filename = [animals{n} '_Rec_avg' suffix{n}];
    cd(['G:\Data\Reach\' animals{n}]);

    load(filename)

    locoC = [];
    reachC = [];
    locoM = [];
    reachM = [];
    for i = 1:numFiles
        locoM = [locoM avg(i).meanMuscles.loco];
        reachM = [reachM avg(i).meanMuscles.reachMvmt]; 
        
        numCells = size(avg(i).meanCells.loco,2);
        for j = 1:numCells
            if mean(avg(i).meanCells.loco(:,j)) > fr_th
                locoC = [locoC avg(i).meanCells.loco(:,j)];
            end
            
            if mean(avg(i).meanCells.reachMvmt(:,j)) > fr_th
                reachC = [reachC avg(i).meanCells.reachMvmt(:,j)];
            end
        end
    end
    locoC = bsxfun(@minus, locoC, mean(locoC));
    reachC = bsxfun(@minus, reachC, mean(reachC));
    locoM = bsxfun(@minus, locoM, mean(locoM));
    reachM = bsxfun(@minus, reachM, mean(reachM));
    
    %First create low-D representations of the data
    [PCsC,modesC,evaluesC] = princomp(reachC);
    varCapReachC = cumsum(evaluesC)./sum(evaluesC);
    varCapReachC(1:PCAdim)

    [PCsM,modesM,evaluesM] = princomp(reachM);
    varCapReachM = cumsum(evaluesM)./sum(evaluesM);
    varCapReachM(1:PCAdim)

%     modesC = bsxfun(@minus, modesC, mean(modesC));
%     modesM = bsxfun(@minus, modesM, mean(modesM));
    
    [A,B,r,U,V] = canoncorr(modesC(:,1:PCAdim),modesM(:,1:PCAdim));
    CorrsReach(:,n) = r;

    %calc frac of OVERALL neural variance captured by each set of canonical coefficients
    %The canonical variables (~modes) are uncorrelated and each have equal variance (for some reason)
    %But their scaling is arbitrary - they are intended to fit the equivalent variable from the other set
    %So we don't calculate the variance capture directly from diag(cov(U)) or diag(cov(V))
    %We also note here that these coefficient vectors aren't quite orthogonal like PCs are
    %We will orthonormalize the coefficient vectors, but in a modified way so that fraction of mode variance captured is preserved in their norm
    %Then we multiply these my the variance across modes. The coefficients are just weights on those modes, so these weight reflect the variance captured too
    %See pic in phone if confused
    
%     onA = Gram_Schmidt_Process_mod(A);
%     onB = Gram_Schmidt_Process_mod(B);
%     
%     NeuralVarCap_Reach_CVs(:,n) = evaluesC(1:PCAdim)'*abs(onA) / sum(evaluesC(1:PCAdim)'*abs(onA))*varCapReachC(PCAdim);
%     MuscVarCap_Reach_CVs(:,n) = evaluesM(1:PCAdim)'*abs(onB) / sum(evaluesM(1:PCAdim)'*abs(onB))*varCapReachM(PCAdim);
    
    onA = Gram_Schmidt_Process(A);
    onB = Gram_Schmidt_Process(B);

%     varFrac = zeros(1,PCAdim);
%     for i = 1:PCAdim
%         CCAmodeC = modesC(:,1:PCAdim)*onA(;
    CCAmodesC = modesC(:,1:PCAdim)*onA;
    CCAmodesM = modesM(:,1:PCAdim)*onB;

%     NeuralVarCap_Reach_CVs(:,n) = (diag(cov(CCAmodesC)) / trace(cov(CCAmodesC))).*(evaluesC(1:PCAdim)/sum(evaluesC));
%     MuscVarCap_Reach_CVs(:,n) = (diag(cov(CCAmodesM)) / trace(cov(CCAmodesM))).*(evaluesM(1:PCAdim)/sum(evaluesM));
%     a = trace(cov(modesC(:,1:PCAdim)))
%     b = trace(cov(CCAmodesC))
%     pause
%     
    NeuralVarCap_Reach_CVs(:,n) = diag(cov(CCAmodesC)) / trace(cov(modesC))
    MuscVarCap_Reach_CVs(:,n) = diag(cov(CCAmodesM)) / trace(cov(modesM))
%     pause
    
%     a = (diag(cov(CCAmodesC)) / trace(cov(CCAmodesC)))
%     b = (evaluesC(1:PCAdim)/sum(evaluesC(1:PCAdim)))
%     a+b
%     pause

        figure(n)
    for i = 1:plotDim
        subplot(1,plotDim,i)
        plot(U(:,i),'r');
        hold on
        plot(V(:,i),'k');
        title(num2str(r(i)));
        xlim([0 450])
        set(gca,'box','off')
        set(gca,'TickDir','out')
        axis square
    end

    %First create low-D representations of the data
    [PCsC,modesC,evaluesC] = princomp(locoC);
    varCapLocoC = cumsum(evaluesC)./sum(evaluesC);
    varCapLocoC(1:PCAdim)

    [PCsM,modesM,evaluesM] = princomp(locoM);
    varCapLocoM = cumsum(evaluesM)./sum(evaluesM);
    varCapLocoM(1:PCAdim)
% 
%     modesC = bsxfun(@minus, modesC, mean(modesC));
%     modesM = bsxfun(@minus, modesM, mean(modesM));
    
    [A,B,r,U,V] = canoncorr(modesC(:,1:PCAdim),modesM(:,1:PCAdim));
    CorrsLoco(:,n) = r;
    
%     onA = Gram_Schmidt_Process_mod(A);
%     onB = Gram_Schmidt_Process_mod(B);
% 
%     NeuralVarCap_Loco_CVs(:,n) = evaluesC(1:PCAdim)'*abs(onA) / sum(evaluesC(1:PCAdim)'*abs(onA))*varCapLocoC(PCAdim);
%     MuscVarCap_Loco_CVs(:,n) = evaluesM(1:PCAdim)'*abs(onB) / sum(evaluesM(1:PCAdim)'*abs(onB))*varCapLocoM(PCAdim);

    onA = Gram_Schmidt_Process(A);
    onB = Gram_Schmidt_Process(B);

    CCAmodesC = modesC(:,1:PCAdim)*onA;
    CCAmodesM = modesM(:,1:PCAdim)*onB;

    NeuralVarCap_Loco_CVs(:,n) = diag(cov(CCAmodesC)) / trace(cov(modesC))
    MuscVarCap_Loco_CVs(:,n) = diag(cov(CCAmodesM)) / trace(cov(modesM))
%     pause

    figure(n+10)
    for i = 1:plotDim
        subplot(1,plotDim,i)
        plot(U(:,i),'r');
        hold on
        plot(V(:,i),'k');
        title(num2str(r(i)));
        xlim([0 450])
        set(gca,'box','off')
        set(gca,'TickDir','out')
        axis square
    end

end

figure(100)
errorbar((1:plotDim)+0.13,mean(CorrsLoco,2),std(CorrsLoco,[],2)/sqrt(length(animals)),'ko-')
hold on
errorbar((1:plotDim)-0.13,mean(CorrsReach,2),std(CorrsReach,[],2)/sqrt(length(animals)),'ro-')
ylim([0 1.25])
xlim([0.5 6.5])
set(gca,'box','off')
set(gca,'TickDir','out')
axis square

cumNeuralVarCap_Reach_CVs = cumsum(NeuralVarCap_Reach_CVs,1);
cumMuscVarCap_Reach_CVs = cumsum(MuscVarCap_Reach_CVs,1);
cumNeuralVarCap_Loco_CVs = cumsum(NeuralVarCap_Loco_CVs,1);
cumMuscVarCap_Loco_CVs = cumsum(MuscVarCap_Loco_CVs,1);

figure(101)
errorbar((1:plotDim)-0.13,mean(cumNeuralVarCap_Reach_CVs,2),std(cumNeuralVarCap_Reach_CVs,[],2)/sqrt(length(animals)),'r-')
hold on
errorbar((1:plotDim)-0.13,mean(cumMuscVarCap_Reach_CVs,2),std(cumMuscVarCap_Reach_CVs,[],2)/sqrt(length(animals)),'k-')
errorbar((1:plotDim)+0.13,mean(cumNeuralVarCap_Loco_CVs,2),std(cumNeuralVarCap_Loco_CVs,[],2)/sqrt(length(animals)),'m-')
errorbar((1:plotDim)+0.13,mean(cumMuscVarCap_Loco_CVs,2),std(cumMuscVarCap_Loco_CVs,[],2)/sqrt(length(animals)),'b-')
xlim([0.5 6.5])
ylim([0 1])
set(gca,'box','off')
set(gca,'TickDir','out')
axis square

mean(sum([NeuralVarCap_Reach_CVs(5:6,:) NeuralVarCap_Loco_CVs(5:6,:)]))
std(sum([NeuralVarCap_Reach_CVs(5:6,:) NeuralVarCap_Loco_CVs(5:6,:)]))/sqrt(6)

mean(sum([MuscVarCap_Reach_CVs(5:6,:) MuscVarCap_Loco_CVs(5:6,:)]))
std(sum([MuscVarCap_Reach_CVs(5:6,:) MuscVarCap_Loco_CVs(5:6,:)]))/sqrt(6)
