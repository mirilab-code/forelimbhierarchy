function [c2r,r2c, c2r_coeffs, r2c_coeffs] = get_CCM_pvalues(best_c2r_coeffs,best_r2c_coeffs, c2r_null_coeffs, r2c_null_coeffs)
%GET_CCM_PVALUES Summary of this function goes here
%   Detailed explanation goes here

c2r_coeffs = best_c2r_coeffs(~isnan(best_c2r_coeffs));
r2c_coeffs = best_r2c_coeffs(~isnan(best_r2c_coeffs));

c2r_null_coeffs = c2r_null_coeffs(~isnan(best_c2r_coeffs),:);
r2c_null_coeffs = r2c_null_coeffs(~isnan(best_r2c_coeffs),:);

c2r = zeros(size(c2r_coeffs));
r2c = zeros(size(r2c_coeffs));

for i=1:size(c2r,1)
    num_null = 300 - sum(isnan(c2r_null_coeffs(i,:)));
    c2r(i,:) = sum(c2r_null_coeffs(i,:) > c2r_coeffs(i)) ./ num_null;
    num_null = 300 - sum(isnan(r2c_null_coeffs(i,:)));
    r2c(i,:) = sum(r2c_null_coeffs(i,:) > r2c_coeffs(i)) ./ num_null;
end


%{
f2 = figure();
hold on
histogram(best_r2c_coeffs(r2c<=0.05), 'Normalization', 'probability', 'Facecolor', 'red', ...
    'Facealpha', .5, BinWidth=0.05);
histogram(best_c2r_coeffs(c2r<=0.05), 'Normalization', 'probability', 'Facecolor', 'blue', ...
    'Facealpha', .5, BinWidth=0.05);
title('ccm coeffs red=rfa->cfa, blue=cfa->rfa')
ylabel('fraction of significant neuron pairs')
xlabel('ccm_coeffs')
hold off

f1 = figure();
hold on
histogram(r2c, 'Normalization', 'probability', 'Facecolor', 'red', ...
    'Facealpha', .5, BinWidth=0.05);
histogram(c2r, 'Normalization', 'probability', 'Facecolor', 'blue', ...
    'Facealpha', .5, BinWidth=0.05);
title('red=rfa->cfa, blue=cfa->rfa')
ylabel('fraction of neuron pairs')
xlabel('p values')
hold off

%}


