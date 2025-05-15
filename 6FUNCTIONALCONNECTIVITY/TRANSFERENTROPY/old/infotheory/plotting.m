%%
totalN = size(peakTE,1);
maxTE = max(peakTE(:));
minTE = min(peakTE(:));

figure(1);
imagesc(peakTE);
ax = gca;
ax.YDir = 'normal';

line([nUnits(1),nUnits(1)], [0,totalN], 'Color', 'r');
line( [0,totalN], [nUnits(1),nUnits(1)], 'Color', 'r');
caxis([minTE maxTE]);

xlabel('CFA | RFA');
ylabel('CFA | RFA');

% for the movement only matrix
maxTE_movement = max(peakTE_movement(:));
minTE_movement = min(peakTE_movement(:));

figure(2);
imagesc(peakTE_movement);
ax = gca;
ax.YDir = 'normal';

line([nUnits(1),nUnits(1)], [0,totalN], 'Color', 'r');
line( [0,totalN], [nUnits(1),nUnits(1)], 'Color', 'r');
caxis([minTE_movement maxTE_movement]);

xlabel('CFA | RFA');
ylabel('CFA | RFA');

%% compare movement vs non movement
difference = peakTE_movement - peakTE;

figure(13);
imagesc(difference);
ax = gca;
ax.YDir = 'normal';

line([nUnits(1),nUnits(1)], [0,totalN], 'Color', 'r');
line( [0,totalN], [nUnits(1),nUnits(1)], 'Color', 'r');
caxis([min(difference(:)) max(difference(:))]);

xlabel('CFA | RFA');
ylabel('CFA | RFA');

dhist = difference(:);
dhist = dhist(dhist~=0);
figure(1);
histogram(dhist);

differenceCR = CFAtoRFA_movement-CFAtoRFA;
differenceRC = RFAtoCFA_movement-RFAtoCFA;

dCRhist = differenceCR(:);
dCRhist = dCRhist(dCRhist~=0);
dRChist = differenceRC(:);
dRChist = dRChist(dRChist~=0);

figure(2);
histogram(dCRhist);
figure(3);
histogram(dRChist);
%%
TEhist = peakTE(:);
TEhist(TEhist==0) = [];
histogram(TEhist)

%% make a matrix of all the delays and plot them
delays = reshape(TEdelays,[],size(TEdelays,3),1);
figure(234);
hold on
for i=1:size(delays,1)
    plot(delays(i,:));
end
hold off

%%
Phist = Pmap(:);
Phist = Phist(~isnan(Phist));

figure(1);
imagesc(Pmap);
ax = gca;
ax.YDir = 'normal';

figure(2);
histogram(Phist, 'BinWidth',0.02);


%%

%%
sigRC = plot_sig(RFAtoCFA_p,33);
RCphist = RFAtoCFA_p(:);
RCphist = RCphist(~isnan(RCphist));

figure(321);
histogram(RCphist);

sigCR = plot_sig(CFAtoRFA_p,34);
CRphist = CFAtoRFA_p(:);
CRphist = CRphist(~isnan(CRphist));

figure(322);
histogram(CRphist);

%%
figure(1);
imagesc(RFAtoCFA_te);
figure(2);
imagesc(CFAtoRFA_te');
figure(3);
imagesc(RFAtoCFA_te - CFAtoRFA_te');





%%
Matrix3DMovie(TEdelays)
