%%



[coef1,score1,latent1] = pca(RFA');
[coef2,score2,latent2] = pca(CFA');

figure;
subplot(2,1,1);
imagesc(coef1'*coef1);

subplot(2,1,2);
imagesc(coef2'*coef2);

%%
nPCs = 25;

score1()

%%

m=18;
n=8;
L1=zeros(m,n);
for j=1:n
    L1(:,j)=((m-n+j-1).*(m-n+j)).^(-1/2).*[ones(m-n+j-1,1) ; -(m-n+j-1) ; zeros(n-j,1)];
end

m=18;
n=5;
L2=zeros(m,n);
for j=1:n
    L2(:,j)=((m-n+j-1).*(m-n+j)).^(-1/2).*[ones(m-n+j-1,1) ; -(m-n+j-1) ; zeros(n-j,1)];
end
%%




%%