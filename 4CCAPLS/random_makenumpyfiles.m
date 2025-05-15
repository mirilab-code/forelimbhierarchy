%%

R = {};
vc_RFA = {};
vc_CFA = {};
for i=1:length(lags)
    lag = lags(i);
    [CV1,CV2,r,score1,score2,varcap1,varcap2] = pca_then_cca(RFA_RG,CFA_RG(:,:,i),0.95);
    R{i} = r;
    vc_RFA{i} = varcap1;
    vc_CFA{i} = varcap2;
end

%%
% cmap = [parula(5); [1 0 0]; flipud(parula(5))];
cmap = jet(length(lags));

figure;
hold on
for i=1:length(lags)
    line = plot(R{i},'Color', cmap(i,:),'DisplayName',num2str(lags(i)));
    row1 = dataTipTextRow('lag',line.XData*0+lags(i));
    line.DataTipTemplate.DataTipRows(end+1) = row1;


end
% legend;




%%


%%
A = D(1).RFA_RG;
B = D(1).CFA_RG;
C = D(7).CFA_RG;

sameR = cca_all_lags(A,B,lags);
diffR = cca_all_lags(A,C,lags);




%%












%%
[X,Y,r,U,V] = canoncorr(CFA_RG(:,:,61)',RFA_RG');

%%
figure;
hold on
plot(r)


figure;
subplot(2,1,1)
imagesc(CFA_RG(:,:,61))
subplot(2,1,2)
imagesc(RFA_RG)

figure;
subplot(2,1,1)
imagesc(CFA_RG(:,:,121))
subplot(2,1,2)
imagesc(RFA_RG)


%%
R = [];
for i=1:size(CFA_RG,3)
    [~,~,r,~,~] = canoncorr(CFA_RG(:,:,i)',RFA_RG');
    R = [R; r];

end
figure
hold on
plot(R', 'Color', [0.4 0.4 0.4 0.7])


%%
figure;
subplot(3,3,1)
imagesc(CFA_RG(:,:,1))
title('CFA lag -60')
subplot(3,3,2)
imagesc(RFA_RG)
title('RFA')
subplot(3,3,3)
plot(R(1,:))

subplot(3,3,4)
imagesc(CFA_RG(:,:,61))
title('CFA lag 0')
subplot(3,3,5)
imagesc(RFA_RG)
title('RFA')
subplot(3,3,6)
plot(R(61,:))

subplot(3,3,7)
imagesc(CFA_RG(:,:,121))
title('CFA lag -60')
subplot(3,3,8)
imagesc(RFA_RG)
title('RFA')
subplot(3,3,9)
plot(R(121,:))

%%
[coef_RFA,score_RFA,latent_RFA] = pca(RFA_RG');
[coef_CFA_0,score_CFA_0,latent_CFA_0] = pca(CFA_RG(:,:,61)');
[coef_CFA_minus60,score_CFA_minus60,latent_CFA_minus60] = pca(CFA_RG(:,:,1)');
%%
v_RFA = cumsum(latent_RFA)/sum(latent_RFA);
v_CFA_0 = cumsum(latent_CFA_0)/sum(latent_CFA_0);
v_CFA_minus60 = cumsum(latent_CFA_minus60)/sum(latent_CFA_minus60);

dim_RFA = find(v_RFA>=0.95,1);
dim_CFA_0 = find(v_CFA_0>=0.95,1);
dim_CFA_minus60 = find(v_CFA_minus60>=0.95,1);

dim = max([dim_RFA dim_CFA_0 dim_CFA_minus60]);

M_RFA = score_RFA(:,1:dim);
M_CFA_0 = score_CFA_0(:,1:dim);
M_CFA_minus60 = score_CFA_minus60(:,1:dim);

%%
[~,~,r1,U1,V1] = canoncorr(M_RFA,M_CFA_0);
[~,~,r2,U2,V2] = canoncorr(M_RFA,M_CFA_minus60);

figure
subplot(2,3,1);
imagesc(M_CFA_0');
title('CFA lag 0')
subplot(2,3,2);
imagesc(M_RFA');
subplot(2,3,3);
plot(r1)

subplot(2,3,4);
imagesc(M_CFA_minus60');
title('CFA lag -60')
subplot(2,3,5);
imagesc(M_RFA');
subplot(2,3,6);
plot(r2)


%%

















%%