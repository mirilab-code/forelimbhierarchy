%%

allunits = load('C:\Users\mirilab\Documents\Mark\CCAPLS\all_CFA_RFA_units.mat');
allunits = allunits.all_CFA_RFA_units;
units_cfa = [];
units_rfa = [];

for i=1:length(allunits)
    this = allunits{i};
    units_cfa = [units_cfa; length(this{1})];
    units_rfa = [units_rfa; length(this{2})];

end

%%
minCFA = min(units_cfa);
minRFA = min(units_rfa);

maxCFA = max(units_cfa);
maxRFA = max(units_rfa);

medCFA = median(units_cfa);
medRFA = median(units_rfa);

fprintf('CFA: %d-%d, median=%d \n',minCFA,maxCFA,medCFA);
fprintf('RFA: %d-%d, median=%d \n',minRFA,maxRFA,medRFA);