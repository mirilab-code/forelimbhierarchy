cfa_inh = cfa_width<13;
cfa_exc = ~cfa_inh;
rfa_inh = rfa_width<13;
rfa_exc = ~rfa_inh;

cfa_fr_exc = cfa_fr(cfa_exc);
rfa_fr_exc = rfa_fr(rfa_exc);
cfa_fr_inh = cfa_fr(cfa_inh);
rfa_fr_inh = rfa_fr(rfa_inh);

%%
cfa_units = 1:length(cfa_fr);
rfa_units = 1:length(rfa_fr);



%%
figure;
hold on
yline(13)
s = scatter(cfa_fr,cfa_width);
row = dataTipTextRow('unit',cfa_units);
s.DataTipTemplate.DataTipRows(end+1) = row;
title('cfa')
xlabel('firing rate');
ylabel('width')

figure;
hold on
yline(13)
s = scatter(rfa_fr,rfa_width);
row = dataTipTextRow('unit',rfa_units);
s.DataTipTemplate.DataTipRows(end+1) = row;
title('rfa')
xlabel('firing rate');
ylabel('width')