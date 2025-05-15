%% analyze 10k pairs
% these are the master variables that will persist throughout loading
% different experiments

RtoC_TE = [];
RtoC_P = [];
CtoR_TE = [];
CtoR_P = [];


%%
NULLTE = cat(3,nullTE,nullTE2);
ALLOWED = cat(3,allowed_shifts,allowed_shifts2);
%
thisP = get_pvals(observedTE,NULLTE,ALLOWED);
sig = thisP<0.05;
sigTE = observedTE;
sigTE(~sig) = nan;

% imagesc(thisP)


% pairs: is the rows and columns to find the pairs on P. not necessarily
% significant but the highest firing rate

RFAind = pairs(:,1);
CFAind = pairs(:,2)+nRFA;

for i=1:length(CFAind)
    n1 = RFAind(i);
    n2 = CFAind(i);
    disp([i n1 n2]);
    RtoCp = thisP(n1,n2);
    RtoC_P = [RtoC_P; RtoCp];
    CtoRp = thisP(n2,n1);
    CtoR_P = [CtoR_P; CtoRp];
    
    RtoCte = observedTE(n1,n2);
    RtoC_TE = [RtoC_TE; RtoCte];
    CtoRte = observedTE(n2,n1);
    CtoR_TE = [CtoR_TE; CtoRte];

end





fprintf('done with %d !\n', length(CFAind));

%% plotting
figure;
tiledlayout(2,1);
nexttile;
hold on
histogram(RtoC_P, 'BinWidth',0.05);
histogram(CtoR_P, 'BinWidth',0.05);
hold off
legend({'RFA to CFA','CFA to RFA'});

nexttile;
hold on
histogram(RtoC_TE, 'BinWidth',0.00001);
histogram(CtoR_TE, 'BinWidth',0.00001);
hold off
legend({'RFA to CFA','CFA to RFA'});





%%
save('TEresults.mat','RtoC_P','CtoR_P','RtoC_TE','CtoR_TE');
