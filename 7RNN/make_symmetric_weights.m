%% take the trained interregional weightsand make a new symmetrical distribution

%% load the connectivity matrix and do stuff with it

addpath(genpath('Z:/code'));
% W = readNPY('newW.npy');
files = dir('newW_*.npy');
Ws = {files.name};

W = zeros(1000,1000,length(Ws));
CtoR = [];
RtoC = [];
for i=1:length(Ws)
    w_file = Ws{i};
    w = readNPY(w_file);
    W(:,:,i) = w;
    cfa_to_rfa = w(501:end,1:500);
    rfa_to_cfa = w(1:500,501:end);
    CtoR = [CtoR; cfa_to_rfa(:)];
    RtoC = [RtoC; rfa_to_cfa(:)];
end
CtoR(CtoR==0) = [];
RtoC(RtoC==0) = [];

%%
oneW = W(:,:,1);

w_intra = oneW;
w_intra(501:end,1:500) = 0;
w_intra(1:500,501:end) = 0;
w_sym_mask = oneW;
w_sym_mask(1:500,1:500) = 0;
w_sym_mask(501:end,501:end) = 0;

Wsym = w_intra;
%% now make the w_sym_mask symmetric
weight_pool = cfa_to_rfa(:);
weight_pool(weight_pool==0) = [];

weight_pool = [weight_pool(randperm(length(weight_pool))) weight_pool(randperm(length(weight_pool)))];
nonz = 0;
for i=1:length(Wsym)
    for j=1:length(Wsym)
        s = w_sym_mask(i,j);
        if(s > 0)
            nonz = nonz+1;
            Wsym(i,j) = weight_pool(nonz);
        end
    end
end


% chec Wsym for symmetry
newRtoC = Wsym(1:500,501:end);
newCtoR = Wsym(501:end,1:500);

newRtoC(newRtoC==0) = [];
newCtoR(newCtoR==0) = [];

figure;
clf("reset")
hold on
boxplot([newRtoC; newCtoR]','orientation', 'horizontal', 'Whisker',Inf);
hold off
annotation('textbox',[.2 .27 .3 .3],'String','cfa->rfa','FitBoxToText','on');
annotation('textbox',[.2 .17 .3 .3],'String','rfa->cfa','FitBoxToText','on');


%%
Wsym_mask = Wsym;
Wsym_mask(Wsym>0) = 1;
Wsym_mask(Wsym<0) = -1;
writeNPY(Wsym,'data/w_sym.npy');
writeNPY(Wsym_mask,'data/w_sym_mask.npy');

