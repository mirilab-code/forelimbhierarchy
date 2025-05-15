% %%%Analyze firing rate data from experiments with both loocomotor and precision reach data
close all
% clear all
format longg

%%% Params
rawRate = 40000;
dsRate = 1000;
subsamp = rawRate/dsRate;
numMuscles = 6;
smoWin = 10; % for std filtering
smoWinFilt = 200;
flexThresh = 0.05;

lp = 40;
smoFiltSTD = 20;

preSampsJoy = 0.4*dsRate;
postSampsJoy = 0.1*dsRate;
preSampsRepos = 0.125*dsRate;
postSampsRepos = 0.375*dsRate;
preSampsMvmt = 0.15*dsRate;
postSampsMvmt = 0.35*dsRate;
    
baselineWidth = 100; %in s
cliqueScale = 1;
numFiles = 11;
cliqueSizes = 25;
check_win = 5; %for outlier removal
check_std_scale_start = 5;
size_for_clique_search = 100;

animal = 'C32';
suffix = '_20ms_filt';

cd(['G:\Data\Reach\' animal]);
load('C32_Rec_20ms')

for i = 1:numFiles
    filename(i).mat = [animal '_Rec_' num2str(i) 'r' char(suffix)];
end

%%Find behavioral events for behavioral triggered averages
for i = 1:numFiles
    i
    
    %Find the joystick repositionings - works
    joystick = data(i).behReach(:,2);
    joystickSmo = smooth(diff(joystick),10,'sgolay',3);
    moveThresh = 0.02;
    
    a =   joystickSmo > moveThresh;
    b = [diff(a); 0];
    ind = 1;
    reposSamps = [];
    while ind <= length(joystickSmo)
        if b(ind) == 1
            reposSamps = [reposSamps; ind];
            ind = ind+5000; %skip 5 s
        else ind = ind+1;
        end
    end
    
    %get rid of problematic edge events
    if reposSamps(1) <  postSampsRepos
        reposSamps(1) = [];
    end
    if reposSamps(end) > data(i).reachSamps - dsRate;
        reposSamps(end) = [];
    end
    
    %Find reposSamps with pulls after
    preWin = 0.200;
    postWin = 1;
    moveThresh = -3*10^-3;
    keep = [];
    startSampsFixed = [];
    for j = 1:length(reposSamps)
        a = joystickSmo(reposSamps(j)+preWin*dsRate:reposSamps(j)+postWin*dsRate);
        b = find(a<moveThresh);
        if ~isempty(b)
            keep = [keep; j];
            startSampsFixed = [startSampsFixed; reposSamps(j) + preWin*dsRate + b(1)];
        end
    end
    reposSamps = reposSamps(keep);
    numJoy = length(startSampsFixed)
    numRepos = length(reposSamps)
    
    %Find muscle threshold crossings
    smoBiSig = smooth(abs(data(i).behReach(:,6)),10,'sgolay',3);
    smoEDCSig = smooth(abs(data(i).behReach(:,8)),10,'sgolay',3);

    smoFlexorSig = smoBiSig+smoEDCSig;
    [x,n] = hist(smoFlexorSig,500);
    [m,ind] = max(x);
    baseline = n(ind);
    smoFlexorSig = smoFlexorSig - baseline;

    preWin = 0.0;
    postWin = 0.5;
    mvmtSamps = [];
    for j = 1:length(reposSamps)
        a = smoFlexorSig(reposSamps(j)+preWin*dsRate:reposSamps(j)+postWin*dsRate);
        b = find(a > flexThresh);
        if ~isempty(b)
            mvmtSamps = [mvmtSamps; reposSamps(j) + preWin*dsRate + b(1)];
        end
    end
    numMvmt = length(mvmtSamps)
     
    %sanity check
    figure(60+i)
    plot(joystickSmo,'r')
    hold on
    plot(data(i).behReach(:,1),'k')
    for j = 1:length(reposSamps)
        plot(reposSamps(j),joystickSmo(reposSamps(j)),'go')
        plot(startSampsFixed(j),joystickSmo(startSampsFixed(j)),'mo')
    end

    %sanity check
    figure(80+i)    
    plot(smoFlexorSig,'r')
    hold on
    plot(data(i).behReach(:,1),'k')
    for j = 1:length(mvmtSamps)
        plot(mvmtSamps(j),smoFlexorSig(mvmtSamps(j)),'go')
    end
    
  
    %compute signalModes from reach data
    signalModes = zeros(1,numMuscles);
    EMG = data(i).behReach(:,4:9);
    for j = 1:numMuscles
       dataBaseline = EMG(round(data(i).reachSamps/2)-round(baselineWidth/2*dsRate):round(data(i).reachSamps/2)+round(baselineWidth/2*dsRate),j);
       dataBaselineSmo = 0*dataBaseline;
       for n = smoWin+1:length(dataBaselineSmo)-smoWin
           dataBaselineSmo(n) = std(dataBaseline(n-smoWin:n+smoWin));
       end
       [x,n] = hist(dataBaselineSmo,500);
       [m,ind] = max(x);
       signalModes(j) = n(ind);
    end
    
    %First for joystick pulls
    %First make smoothed EMG arrays FROM RAW UNRECTIFIED DATA
    joys = zeros([numJoy preSampsJoy+postSampsJoy+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numJoy
           raw = EMG(startSampsFixed(k)-preSampsJoy-smoWinFilt:startSampsFixed(k)+postSampsJoy+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           joys(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers to reduce number of trials to search (speed up clique finding)
    check_std_scale = check_std_scale_start;
    while numJoy > size_for_clique_search
        mean_joys = squeeze(mean(joys));
        std_joys = squeeze(std(joys,0,1));
        elim = [];
        for p = 1:numJoy
            for j = check_win+1:size(joys,2)-check_win
                for k = 1:numMuscles
                    if abs(mean(joys(p,j-check_win:j+check_win,k)) - mean(mean_joys(j-check_win:j+check_win,k))) > check_std_scale*mean(std_joys(j-check_win:j+check_win,k));
                        elim = [elim; p];
                    end
                end
            end
        end
        elim = unique(elim);
        joys(elim,:,:) = [];
        numJoy = size(joys,1)
        if numJoy > 1.2*size_for_clique_search
            check_std_scale = check_std_scale-0.5
        else check_std_scale = check_std_scale-0.25
        end
    end
           
    %Now calculate pairwise distances
    distsAll = zeros(numJoy,numJoy,numMuscles);
    for j = 1:numMuscles
       distsAll(:,:,j) = squareform(pdist(joys(:,:,j),'correlation')*cliqueScale); %to make them more than 100
    end

    %Now collapse to 2D array
    dists = zeros(numJoy,numJoy);
    for j = 1:numJoy
       for k = 1:numJoy
           dist = 0;
           for m = 1:numMuscles
               dist = dist + distsAll(j,k,m)^2;
           end
           dist = sqrt(dist);
           dists(j,k) = dist;
       end
    end
    
    %now solve the clique problem for different sizes
    distThreshold = 1;
    maxCliqueSize = 0;
    for c = 1:length(cliqueSizes)
        cliqueSizeCurrent = cliqueSizes(c);
        while maxCliqueSize < cliqueSizeCurrent
           distsBin = dists < distThreshold;
           for j = 1:numJoy
               distsBin(j,j) = 0;
           end
           MC = maximalCliques(distsBin,'v2');
           sizes = sum(MC);
           [maxCliqueSize,ind] = max(sizes);
           if maxCliqueSize > cliqueSizeCurrent/2
               maxCliqueSize
           end
           distThreshold = distThreshold*1.03;
        end
        goodJoys = joys(find(MC(:,ind)),:,:);
        goodJoySamps = startSampsFixed(find(MC(:,ind)));
        numGoodJoys = length(find(MC(:,ind)));

        %Make mean muscle activities
        meanMuscles(c).reachJoy = squeeze(mean(goodJoys));
        trialsMuscles(c).reachJoy = goodJoys;
        for j = 1:numMuscles
            figure(i); subplot(6,1,j); plot(meanMuscles(c).reachJoy(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsJoy+postSampsJoy+1 length(data(i).units)]);
        trialsCells(c).reachJoy = zeros([numGoodJoys preSampsJoy+postSampsJoy+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsJoy+postSampsJoy+1 numGoodJoys]);
            for k = 1:numGoodJoys
                mat(:,k) = data(i).units(j).frReach(goodJoySamps(k)-preSampsJoy:goodJoySamps(k)+postSampsJoy);
            end
            trialsCells(c).reachJoy(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachJoy = meanC;
    end
    
    %Now for Repositionings
    %First make smoothed EMG arrays FROM RAW UNRECTIFIED DATA
    Repos = zeros([numRepos preSampsRepos+postSampsRepos+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numRepos
           raw = EMG(reposSamps(k)-preSampsRepos-smoWinFilt:reposSamps(k)+postSampsRepos+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           Repos(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers to reduce number of trials to search (speed up clique finding)
    check_std_scale = check_std_scale_start;
    while numRepos > size_for_clique_search
        mean_Repos = squeeze(mean(Repos));
        std_Repos = squeeze(std(Repos,0,1));
        elim = [];
        for p = 1:numRepos
            for j = check_win+1:size(Repos,2)-check_win
                for k = 1:numMuscles
                    if abs(mean(Repos(p,j-check_win:j+check_win,k)) - mean(mean_Repos(j-check_win:j+check_win,k))) > check_std_scale*mean(std_Repos(j-check_win:j+check_win,k));
                        elim = [elim; p];
                    end
                end
            end
        end
        elim = unique(elim);
        Repos(elim,:,:) = [];
        numRepos = size(Repos,1)
        if numRepos > 1.2*size_for_clique_search
            check_std_scale = check_std_scale-0.5
        else check_std_scale = check_std_scale-0.25
        end
    end
    
    %Now calculate pairwise distances
    distsAll = zeros(numRepos,numRepos,numMuscles);
    for j = 1:numMuscles
       distsAll(:,:,j) = squareform(pdist(Repos(:,:,j),'correlation')*cliqueScale); %to make them more than 100
    end

    %Now collapse to 2D array
    dists = zeros(numRepos,numRepos);
    for j = 1:numRepos
       for k = 1:numRepos
           dist = 0;
           for m = 1:numMuscles
               dist = dist + distsAll(j,k,m)^2;
           end
           dist = sqrt(dist);
           dists(j,k) = dist;
       end
    end
    
    %now solve the clique problem for different sizes
    distThreshold = 1;
    maxCliqueSize = 0;
    for c = 1:length(cliqueSizes)
        cliqueSizeCurrent = cliqueSizes(c);
        while maxCliqueSize < cliqueSizeCurrent
           distsBin = dists < distThreshold;
           for j = 1:numRepos
               distsBin(j,j) = 0;
           end
           MC = maximalCliques(distsBin,'v2');
           sizes = sum(MC);
           [maxCliqueSize,ind] = max(sizes);
           if maxCliqueSize > cliqueSizeCurrent/2
               maxCliqueSize
           end
           distThreshold = distThreshold*1.03;
        end
        goodRepos = Repos(find(MC(:,ind)),:,:);
        goodReposSamps = reposSamps(find(MC(:,ind)));
        numGoodRepos = length(find(MC(:,ind)));

        %Make mean muscle activities
        meanMuscles(c).reachRepos = squeeze(mean(goodRepos));
        trialsMuscles(c).reachRepos = goodRepos;
        for j = 1:numMuscles
            figure(20+i); subplot(6,1,j); plot(meanMuscles(c).reachRepos(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsRepos+postSampsRepos+1 length(data(i).units)]);
        trialsCells(c).reachRepos = zeros([numGoodRepos preSampsRepos+postSampsRepos+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsRepos+postSampsRepos+1 numGoodRepos]);
            for k = 1:numGoodRepos
                mat(:,k) = data(i).units(j).frReach(goodReposSamps(k)-preSampsRepos:goodReposSamps(k)+postSampsRepos);
            end
            trialsCells(c).reachRepos(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachRepos = meanC;
    end
   
    %Now for Mvmt
    %First make smoothed EMG arrays FROM RAW UNRECTIFIED DATA
    Mvmt = zeros([numMvmt preSampsMvmt+postSampsMvmt+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numMvmt
           raw = EMG(mvmtSamps(k)-preSampsMvmt-smoWinFilt:mvmtSamps(k)+postSampsMvmt+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           Mvmt(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers to reduce number of trials to search (speed up clique finding)
    check_std_scale = check_std_scale_start;
    while numMvmt > size_for_clique_search
        mean_Mvmt = squeeze(mean(Mvmt));
        std_Mvmt = squeeze(std(Mvmt,0,1));
        elim = [];
        for p = 1:numMvmt
            for j = check_win+1:size(Mvmt,2)-check_win
                for k = 1:numMuscles
                    if abs(mean(Mvmt(p,j-check_win:j+check_win,k)) - mean(mean_Mvmt(j-check_win:j+check_win,k))) > check_std_scale*mean(std_Mvmt(j-check_win:j+check_win,k));
                        elim = [elim; p];
                    end
                end
            end
        end
        elim = unique(elim);
        Mvmt(elim,:,:) = [];
        numMvmt = size(Mvmt,1)
        if numMvmt > 1.2*size_for_clique_search
            check_std_scale = check_std_scale-0.5
        else check_std_scale = check_std_scale-0.25
        end
    end
    
    %Now calculate pairwise distances
    distsAll = zeros(numMvmt,numMvmt,numMuscles);
    for j = 1:numMuscles
       distsAll(:,:,j) = squareform(pdist(Mvmt(:,:,j),'correlation')*cliqueScale); %to make them more than 100
    end

    %Now collapse to 2D array
    dists = zeros(numMvmt,numMvmt);
    for j = 1:numMvmt
       for k = 1:numMvmt
           dist = 0;
           for m = 1:numMuscles
               dist = dist + distsAll(j,k,m)^2;
           end
           dist = sqrt(dist);
           dists(j,k) = dist;
       end
    end
    
    %now solve the clique problem for different sizes
    distThreshold = 1;
    maxCliqueSize = 0;
    for c = 1:length(cliqueSizes)
        cliqueSizeCurrent = cliqueSizes(c);
        while maxCliqueSize < cliqueSizeCurrent
           distsBin = dists < distThreshold;
           for j = 1:numMvmt
               distsBin(j,j) = 0;
           end
           MC = maximalCliques(distsBin,'v2');
           sizes = sum(MC);
           [maxCliqueSize,ind] = max(sizes);
           if maxCliqueSize > cliqueSizeCurrent/2
               maxCliqueSize
           end
           distThreshold = distThreshold*1.03;
        end
        goodMvmt = Mvmt(find(MC(:,ind)),:,:);
        goodMvmtSamps = mvmtSamps(find(MC(:,ind)));
        numGoodMvmt = length(find(MC(:,ind)));

        %Make mean muscle activities
        meanMuscles(c).reachMvmt = squeeze(mean(goodMvmt));
        trialsMuscles(c).reachMvmt = goodMvmt;
        for j = 1:numMuscles
            figure(40+i); subplot(6,1,j); plot(meanMuscles(c).reachMvmt(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsMvmt+postSampsMvmt+1 length(data(i).units)]);
        trialsCells(c).reachMvmt = zeros([numGoodMvmt preSampsMvmt+postSampsMvmt+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsMvmt+postSampsMvmt+1 numGoodMvmt]);
            for k = 1:numGoodMvmt
                mat(:,k) = data(i).units(j).frReach(goodMvmtSamps(k)-preSampsMvmt:goodMvmtSamps(k)+postSampsMvmt);
            end
            trialsCells(c).reachMvmt(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachMvmt = meanC;
    end
    
    sM = size(meanMuscles(c).reachRepos)
    sC = size(meanCells(c).reachRepos)
    save(filename(i).mat,'meanMuscles','meanCells','cliqueSizes')
    disp(['saved ' filename(i).mat])

end

%make grand template
templateJoyMat = zeros(numFiles,size(meanMuscles(end).reachJoy,1),size(meanMuscles(end).reachJoy,2));
templateReposMat = zeros(numFiles,size(meanMuscles(end).reachRepos,1),size(meanMuscles(end).reachRepos,2));
templateMvmtMat = zeros(numFiles,size(meanMuscles(end).reachMvmt,1),size(meanMuscles(end).reachMvmt,2));
for i = 1:numFiles
    load(filename(i).mat,'meanMuscles')
    templateJoyMat(i,:,:) = meanMuscles(end).reachJoy;
    templateReposMat(i,:,:) = meanMuscles(end).reachRepos;
    templateMvmtMat(i,:,:) = meanMuscles(end).reachMvmt;
end
templateJoy = squeeze(mean(templateJoyMat));
templateRepos = squeeze(mean(templateReposMat));
templateMvmt = squeeze(mean(templateMvmtMat));

save(['templates_' animal],'templateJoy','templateRepos','templateMvmt');

for j = 1:numMuscles
    figure(200); 
    
    subplot(6,3,j);
    plot(templateJoy(:,j));
    
    subplot(6,3,j+6);
    plot(templateRepos(:,j));
   
    subplot(6,3,j+12);
    plot(templateMvmt(:,j));

    set(gcf,'Position',[0 0 1200 1200])
end

cliqueSizes = 40;
check_std_scale = 5;
for i = 1:numFiles
    i
    
    %Find the joystick repositionings - works
    joystick = data(i).behReach(:,2);
    joystickSmo = smooth(diff(joystick),10,'sgolay',3);
    moveThresh = 0.02;
    
    a =   joystickSmo > moveThresh;
    b = [diff(a); 0];
    ind = 1;
    reposSamps = [];
    while ind <= length(joystickSmo)
        if b(ind) == 1
            reposSamps = [reposSamps; ind];
            ind = ind+5000; %skip 5 s
        else ind = ind+1;
        end
    end
    
    %get rid of problematic edge events
    if reposSamps(1) <  postSampsRepos
        reposSamps(1) = [];
    end
    if reposSamps(end) > data(i).reachSamps - dsRate;
        reposSamps(end) = [];
    end
    
    %Find reposSamps with pulls after
    preWin = 0.200;
    postWin = 1;
    moveThresh = -3*10^-3;
    keep = [];
    startSampsFixed = [];
    for j = 1:length(reposSamps)
        a = joystickSmo(reposSamps(j)+preWin*dsRate:reposSamps(j)+postWin*dsRate);
        b = find(a<moveThresh);
        if ~isempty(b)
            keep = [keep; j];
            startSampsFixed = [startSampsFixed; reposSamps(j) + preWin*dsRate + b(1)];
        end
    end
    reposSamps = reposSamps(keep);
    numJoy = length(startSampsFixed)
    numRepos = length(reposSamps)
    
    %Find muscle threshold crossings
    smoBiSig = smooth(abs(data(i).behReach(:,6)),10,'sgolay',3);
    smoEDCSig = smooth(abs(data(i).behReach(:,8)),10,'sgolay',3);

    smoFlexorSig = smoBiSig+smoEDCSig;
    [x,n] = hist(smoFlexorSig,500);
    [m,ind] = max(x);
    baseline = n(ind);
    smoFlexorSig = smoFlexorSig - baseline;
    
    preWin = 0.0;
    postWin = 0.5;
    mvmtSamps = [];
    for j = 1:length(reposSamps)
        a = smoFlexorSig(reposSamps(j)+preWin*dsRate:reposSamps(j)+postWin*dsRate);
        b = find(a > flexThresh);
        if ~isempty(b)
            mvmtSamps = [mvmtSamps; reposSamps(j) + preWin*dsRate + b(1)];
        end
    end
    numMvmt = length(mvmtSamps)
    
    %sanity check
    figure(100+i)
    plot(joystickSmo,'r')
    hold on
    plot(data(i).behReach(:,1),'k')
    for j = 1:length(reposSamps)
        plot(reposSamps(j),joystickSmo(reposSamps(j)),'go')
        plot(startSampsFixed(j),joystickSmo(startSampsFixed(j)),'mo')
    end

    %sanity check
    figure(120+i)
    plot(smoFlexorSig,'r')
    hold on
    plot(data(i).behReach(:,1),'k')
    for j = 1:length(mvmtSamps)
        plot(mvmtSamps(j),smoFlexorSig(mvmtSamps(j)),'go')
    end

    %compute signalModes from reach data
    signalModes = zeros(1,numMuscles);
    EMG = data(i).behReach(:,4:9);
    for j = 1:numMuscles
       dataBaseline = EMG(round(data(i).reachSamps/2)-round(baselineWidth/2*dsRate):round(data(i).reachSamps/2)+round(baselineWidth/2*dsRate),j);
       dataBaselineSmo = 0*dataBaseline;
       for n = smoWin+1:length(dataBaselineSmo)-smoWin
           dataBaselineSmo(n) = std(dataBaseline(n-smoWin:n+smoWin));
       end
       [x,n] = hist(dataBaselineSmo,500);
       [m,ind] = max(x);
       signalModes(j) = n(ind);
    end
    
    %First for joystick pulls
    %First make smoothed EMG arrays FROM RAW UNRECTIFIED DATA
    joys = zeros([numJoy preSampsJoy+postSampsJoy+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numJoy
           raw = EMG(startSampsFixed(k)-preSampsJoy-smoWinFilt:startSampsFixed(k)+postSampsJoy+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           joys(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers
    mean_joys = squeeze(mean(joys));
    std_joys = squeeze(std(joys,0,1));
    elim = [];
    for p = 1:numJoy
        for j = check_win+1:size(joys,2)-check_win
            for k = 1:numMuscles
                if abs(mean(joys(p,j-check_win:j+check_win,k)) - mean(mean_joys(j-check_win:j+check_win,k))) > check_std_scale*mean(std_joys(j-check_win:j+check_win,k));
                    elim = [elim; p];
                end
            end
        end
    end
    elim = unique(elim);
    joys(elim,:,:) = [];
    numJoy = size(joys,1)
           
    %Find the traces most similar to template
    dists = zeros(size(joys,1),numMuscles);
    for j = 1:numMuscles
        dists_all = pdist([templateJoy(:,j)'; joys(:,:,j)],'correlation');
        dists(:,j) = dists_all(1:size(joys,1));
    end
    rms_dists =sqrt(sum(dists.^2,2));
    sort_mat = sortrows([rms_dists [1:size(joys,1)]'],1);

    for c = 1:length(cliqueSizes)
        if size(sort_mat,1) > cliqueSizes(c)
            keep = sort_mat(1:cliqueSizes(c),2);
        else keep = sort_mat(:,2);
        end
        goodJoys = joys(keep,:,:);
        goodJoySamps = startSampsFixed(keep);
        numGoodJoys = size(goodJoys,1)

        %Make mean muscle activities
        meanMuscles(c).reachJoy = squeeze(mean(goodJoys));
        trialsMuscles(c).reachJoy = goodJoys;
        for j = 1:numMuscles
            figure(i); subplot(6,1,j); plot(meanMuscles(c).reachJoy(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsJoy+postSampsJoy+1 length(data(i).units)]);
        trialsCells(c).reachJoy = zeros([numGoodJoys preSampsJoy+postSampsJoy+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsJoy+postSampsJoy+1 numGoodJoys]);
            for k = 1:numGoodJoys
                mat(:,k) = data(i).units(j).frReach(goodJoySamps(k)-preSampsJoy:goodJoySamps(k)+postSampsJoy);
            end
            trialsCells(c).reachJoy(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachJoy = meanC;
    end
    
    %Now for Repositionings
    Repos = zeros([numRepos preSampsRepos+postSampsRepos+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numRepos
           raw = EMG(reposSamps(k)-preSampsRepos-smoWinFilt:reposSamps(k)+postSampsRepos+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           Repos(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers
    mean_Repos = squeeze(mean(Repos));
    std_Repos = squeeze(std(Repos,0,1));
    elim = [];
    for p = 1:numRepos
        for j = check_win+1:size(Repos,2)-check_win
            for k = 1:numMuscles
                if abs(mean(Repos(p,j-check_win:j+check_win,k)) - mean(mean_Repos(j-check_win:j+check_win,k))) > check_std_scale*mean(std_Repos(j-check_win:j+check_win,k));
                    elim = [elim; p];
                end
            end
        end
    end
    elim = unique(elim);
    Repos(elim,:,:) = [];
    numRepos = size(Repos,1)
    
    %Find the traces most similar to template
    dists = zeros(size(Repos,1),numMuscles);
    for j = 1:numMuscles
        dists_all = pdist([templateRepos(:,j)'; Repos(:,:,j)],'correlation');
        dists(:,j) = dists_all(1:size(Repos,1));
    end
    rms_dists =sqrt(sum(dists.^2,2));
    sort_mat = sortrows([rms_dists [1:size(Repos,1)]'],1);

    for c = 1:length(cliqueSizes)
        if size(sort_mat,1) > cliqueSizes(c)
            keep = sort_mat(1:cliqueSizes(c),2);
        else keep = sort_mat(:,2);
        end
        goodRepos = Repos(keep,:,:);
        goodReposSamps = reposSamps(keep);
        numGoodRepos = size(goodRepos,1)
        
        %Make mean muscle activities
        meanMuscles(c).reachRepos = squeeze(mean(goodRepos));
        trialsMuscles(c).reachRepos = goodRepos;
        for j = 1:numMuscles
            figure(20+i); subplot(6,1,j); plot(meanMuscles(c).reachRepos(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsRepos+postSampsRepos+1 length(data(i).units)]);
        trialsCells(c).reachRepos = zeros([numGoodRepos preSampsRepos+postSampsRepos+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsRepos+postSampsRepos+1 numGoodRepos]);
            for k = 1:numGoodRepos
                mat(:,k) = data(i).units(j).frReach(goodReposSamps(k)-preSampsRepos:goodReposSamps(k)+postSampsRepos);
            end
            trialsCells(c).reachRepos(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachRepos = meanC;
    end
 
    %Now for Mvmt
    Mvmt = zeros([numMvmt preSampsMvmt+postSampsMvmt+1 numMuscles]);
    for j = 1:numMuscles
       for k = 1:numMvmt
           raw = EMG(mvmtSamps(k)-preSampsMvmt-smoWinFilt:mvmtSamps(k)+postSampsMvmt+smoWinFilt,j);
           rawSmo = filterEMG(raw,lp,smoFiltSTD);
           Mvmt(k,:,j) = rawSmo(smoWinFilt+1:length(raw)-smoWinFilt)-signalModes(j);
       end
    end
    
    %Now remove outliers
    mean_Mvmt = squeeze(mean(Mvmt));
    std_Mvmt = squeeze(std(Mvmt,0,1));
    elim = [];
    for p = 1:numMvmt
        for j = check_win+1:size(Mvmt,2)-check_win
            for k = 1:numMuscles
                if abs(mean(Mvmt(p,j-check_win:j+check_win,k)) - mean(mean_Mvmt(j-check_win:j+check_win,k))) > check_std_scale*mean(std_Mvmt(j-check_win:j+check_win,k));
                    elim = [elim; p];
                end
            end
        end
    end
    elim = unique(elim);
    Mvmt(elim,:,:) = [];
    numMvmt = size(Mvmt,1)
    
    %Find the traces most similar to template
    dists = zeros(size(Mvmt,1),numMuscles);
    for j = 1:numMuscles
        dists_all = pdist([templateMvmt(:,j)'; Mvmt(:,:,j)],'correlation');
        dists(:,j) = dists_all(1:size(Mvmt,1));
    end
    rms_dists =sqrt(sum(dists.^2,2));
    sort_mat = sortrows([rms_dists [1:size(Mvmt,1)]'],1);

    for c = 1:length(cliqueSizes)
        if size(sort_mat,1) > cliqueSizes(c)
            keep = sort_mat(1:cliqueSizes(c),2);
        else keep = sort_mat(:,2);
        end
        goodMvmt = Mvmt(keep,:,:);
        goodMvmtSamps = mvmtSamps(keep);
        numGoodMvmt = size(goodMvmt,1)
        
        %Make mean muscle activities
        meanMuscles(c).reachMvmt = squeeze(mean(goodMvmt));
        trialsMuscles(c).reachMvmt = goodMvmt;
        for j = 1:numMuscles
            figure(40+i); subplot(6,1,j); plot(meanMuscles(c).reachMvmt(:,j));
            title(num2str(i))
            set(gcf,'Position',[0 0 400 1200])
        end

        %Make mean cell activities
        meanC = zeros([preSampsMvmt+postSampsMvmt+1 length(data(i).units)]);
        trialsCells(c).reachMvmt = zeros([numGoodMvmt preSampsMvmt+postSampsMvmt+1 length(data(i).units)]);
        for j = 1:length(data(i).units)
            mat = zeros([preSampsMvmt+postSampsMvmt+1 numGoodMvmt]);
            for k = 1:numGoodMvmt
                mat(:,k) = data(i).units(j).frReach(goodMvmtSamps(k)-preSampsMvmt:goodMvmtSamps(k)+postSampsMvmt);
            end
            trialsCells(c).reachMvmt(:,:,j) = mat';
            meanC(:,j) = mean(mat,2);
        end
        meanCells(c).reachMvmt = meanC;
    end
    
    sM = size(meanMuscles(c(end)).reachRepos)
    sC = size(meanCells(c(end)).reachRepos)
    data_session = data(i);
    save(filename(i).mat,'data_session','meanMuscles','meanCells','trialsMuscles','trialsCells','goodReposSamps','goodJoySamps','goodMvmtSamps','cliqueSizes',...
        'preSampsJoy','postSampsJoy','preSampsRepos','postSampsRepos','preSampsMvmt','postSampsMvmt','reposSamps','startSampsFixed','mvmtSamps')
    disp(['saved ' filename(i).mat])

end

