function [X,Y,r,U,V,varcap1,varcap2] = pca_then_cca(M1,M2,nPCs)

M1 = M1 - mean(M1,2);
M2 = M2 - mean(M2,2);
M1 = M1';
M2 = M2';

dim = nPCs;

[coef1,score1,latent1] = pca(M1);
[coef2,score2,latent2] = pca(M2);

totalv1 = latent1 / sum(latent1);
totalv1 = totalv1(1:dim);
totalv2 = latent2 / sum(latent2);
totalv2 = totalv2(1:dim);

% v1 = latent1(1:dim);
% v2 = latent2(1:dim);
v1 = latent1;
v2 = latent2;
v1 = v1 / sum(v1);
v2 = v2 / sum(v2);
v1 = v1(1:dim);
v2 = v2(1:dim);

X1 = score1;
X2 = score2; 
X1 = score1(:,1:dim);
X2 = score2(:,1:dim);


[X,Y,r,U,V] = canoncorr(X1,X2);

% orthonormal stuff and getting var capture
% orth1 = orth(X);
% Q = X1*orth1;
% varcap1 = flipud(diag(cov(Q)) / trace(cov(X1)));
% 
% orth2 = orth(Y);
% Q = X2*orth2;
% varcap2 = flipud(diag(cov(Q)) / trace(cov(X2)));
cvs_to_pcs1 = X * coef1(:,1:dim)';
orth_cca1 = orth(cvs_to_pcs1');
Q_cca1 = M1 * orth_cca1;
cca_vc1 = flipud(diag(cov(Q_cca1)) / trace(cov(M1)));

cvs_to_pcs2 = Y * coef2(:,1:dim)';
orth_cca2 = orth(cvs_to_pcs2');
Q_cca2 = M2 * orth_cca2;
cca_vc2 = flipud(diag(cov(Q_cca2)) / trace(cov(M2)));

varcap1.pca = v1;
varcap1.cca = cca_vc1;
varcap2.pca = v2;
varcap2.cca = cca_vc2;



end