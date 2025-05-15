%% calculate entropy of each neuron
H = [];
for i=1:length(info)
    m = max(ISI{i});
    pi = probdists{i};
    x = 0:m;
    y = pdf(pi,x);
    y(y==Inf) = 0;
    ii = info{i};
    h = sum(y.*ii);
    
    H = [H; h];    
end
disp('done');


pdparamsH = [distparams H];

%% plot bit rate
figure(421);
cmap = spring(size(distparams,1));
colormap(spring);
g = gscatter(distparams(:,1),distparams(:,2),bits_per_second,cmap);
% set(g,'Marker','square');
xlabel('parameter a (shape)');
ylabel('parameter b (scale)');
cb = colorbar;
ylabel(cb, 'Bits per second');
caxis([min(bits_per_second) max(bits_per_second)]);

leg = legend('example');
set(leg,'visible','off');



%% plot entropy
figure(31);
cmap = spring(size(distparams,1));
colormap(spring);
g = gscatter(distparams(:,1),distparams(:,2),H,cmap);
% set(g,'Marker','square');
xlabel('parameter a (shape)');
ylabel('parameter b (scale)');
cb = colorbar;
ylabel(cb, 'Entropy (bits)');
caxis([min(H) max(H)]);

leg = legend('example');
set(leg,'visible','off');


%% same thing but with firingrates
figure(32);
cmap = spring(size(distparams,1));
colormap(spring);
g = gscatter(distparams(:,1),distparams(:,2),firingrates,cmap);
% set(g,'Marker','square');
xlabel('parameter a (shape)');
ylabel('parameter b (scale)');
cb = colorbar;
ylabel(cb, 'firing rate (Hz)');
caxis([min(firingrates) max(firingrates)]);

leg = legend('example');
set(leg,'visible','off');

%%
figure(33);
tiledlayout(1,2)
nexttile;
scatter(firingrates,H, 'filled');
xlabel('firing rate (Hz)');
ylabel('entropy (bits)');
nexttile;
scatter(firingrates,H, 'filled');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('log(firing rate) (Hz)');
ylabel('log(entropy) (bits)');
%% KL Divergence

KLd = zeros(N,N);
for i=1:N
    disp(i);
    for j=1:N
        m = min(max(ISI{i}),max(ISI{j}));
        P = pdf(probdists{i},1:m);
        Q = pdf(probdists{j},1:m);
        x = sum(P.*log2(P./Q));
        KLd(i,j) = x;
    end
end
imagesc(KLd);


%%


u = u-1;
x = 1:max(ISI{u});
y = pdf(probdists{u},x);

figure(121151);
plot(x,y)
figure(13100);
loglog(x,y)

%%
figure(98)
histogram(running_info{98})
figure(161)
histogram(running_info{161})
figure(126)
histogram(running_info{126})
figure(117)
histogram(running_info{117})


%% make joint prob??????
u1 = 123;
u2 = 65;

ri1 = running_info{u1};
ri2 = running_info{u2};

minlen = min(length(ri1),length(ri2));
ri1 = ri1(1:minlen);
ri2 = ri2(1:minlen);

hist3([ri1 ri2],'BinWidth',1 ,'CDataMode','auto','FaceColor','interp');









%%
