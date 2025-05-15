%%
Pmap = get_pvals(observedTE,nullTE,allowed_shifts);
sig = plot_sig(Pmap,nUnits);

[CFAtoRFA_p,RFAtoCFA_p] = separate(Pmap,nUnits);
[CFAtoRFA_te,RFAtoCFA_te] = separate(observedTE,nUnits);

%% FDR
[h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(Pmap,0.3,'pdep','yes');
figure;
imagesc(h);
ax = gca;
ax.YDir = 'normal';

[CFAtoRFA_q,RFAtoCFA_q] = separate(h,nUnits);


%% save the important numbers!!!!
cd(data_path);
cd ..
if ~exist('TransferEntropy', 'dir')
    mkdir('TransferEntropy')
end
cd('TransferEntropy');


save('TEnumbers.mat','Pmap','CFAtoRFA_p','RFAtoCFA_p','CFAtoRFA_te','RFAtoCFA_te','h','crit_p','adj_ci_cvrg','adj_p','CFAtoRFA_q','RFAtoCFA_q');
disp('done saving!');
%%

































%% everything below this is just me goofing around

% %%
% figure(1);
% histogram(CFAtoRFA_p(:));
% figure(2);
% histogram(RFAtoCFA_p(:));
% 
% 
% %% visualize difference by making the pvals way negative or positive
% CR = CFAtoRFA_q;
% RC = RFAtoCFA_q;
% 
% Qdiff = CR-RC';
% imagesc(Qdiff);
% ax = gca;
% ax.YDir = 'normal';
% 
% Qhist = Qdiff(:);
% Qhist = Qhist(Qhist~=0);
% histogram(Qhist,'BinWidth',0.5);
% 
% 
% %% get convergence/divergence
% 
% [CR_conv,RC_conv,CR_div,RC_div] = convdiv(CFAtoRFA_q,RFAtoCFA_q);
% % CR_conv: how many CFA connections to RFA neurons there are for each RFA neuron, on average
% % RC_conv: how many RFA connections to CFA neurons there are for each CFA neuron, on average
% % CR_div:  how many connections each CFA neuron makes to RFA neurons, on average
% % RC_div:  how many connections each RFA neuron makes to CFA neurons, on average
% 
% % when done with 0s and 1s (using the h matrix which is adjusted
% %   signficance, these numbers are actually percentages
% 
% 
% %%
% figure;
% histogram(RC_div,'BinWidth',0.025);
% title('divergence of RFA to CFA');
% 
% %%
% figure;
% histogram(CR_div,'BinWidth',0.025);
% % plot(CR_div);
% 
% figure;
% histogram(RC_conv,'BinWidth',0.025);
% % plot(RC_conv);
% figure;
% histogram(CR_conv,'BinWidth',0.025);
% % plot(CR_conv);
% 
