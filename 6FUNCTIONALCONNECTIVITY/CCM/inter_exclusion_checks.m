%
mean_FR = zeros(10000,2);
pb = progressbar(' getting firing rate... ');
for i=1:size(top10k_pairs, 1)
    pb.print(i, size(top10k_pairs, 1));
    [CFA_FR, RFA_FR, reaches] = grab_pair(top10k_pairs(i,:));
    [~, bouts] = get_reaching_bouts(reaches, size(CFA_FR,2), -300, 0);
    mean_FR(i, 1) = mean(CFA_FR(bouts));
    mean_FR(i, 2) = mean(RFA_FR(bouts));
end
%save('mean_FR.mat', 'mean_FR')
%
temp_best_c2r_coeffs_keep = best_c2r_coeffs(~isnan(best_c2r_coeffs));
temp_best_r2c_coeffs_keep = best_r2c_coeffs(~isnan(best_c2r_coeffs));
temp_c2r_null_coeffs = c2r_null_coeffs(~isnan(best_c2r_coeffs),:);
temp_r2c_null_coeffs = r2c_null_coeffs(~isnan(best_c2r_coeffs),:);
temp_mean_FR_keep = mean_FR(~isnan(best_c2r_coeffs), :);

best_r2c_coeffs_keep = temp_best_r2c_coeffs_keep(~isnan(temp_best_r2c_coeffs_keep));
best_c2r_coeffs_keep = temp_best_c2r_coeffs_keep(~isnan(temp_best_r2c_coeffs_keep));
c2r_null_coeffs_keep = temp_c2r_null_coeffs(~isnan(temp_best_r2c_coeffs_keep),:);
r2c_null_coeffs_keep = temp_r2c_null_coeffs(~isnan(temp_best_r2c_coeffs_keep),:);
mean_FR_keep = temp_mean_FR_keep(~isnan(temp_best_r2c_coeffs_keep), :);


temp = min(mean_FR_keep,[], 2) > 0;
mean_FR_keep = mean_FR_keep(temp, :);
best_c2r_coeffs_keep = best_c2r_coeffs_keep(temp);
best_r2c_coeffs_keep = best_r2c_coeffs_keep(temp);
c2r_null_coeffs_keep = c2r_null_coeffs_keep(temp, :);
r2c_null_coeffs_keep = r2c_null_coeffs_keep(temp, :);

c2r = zeros(size(best_c2r_coeffs_keep));
r2c = zeros(size(best_r2c_coeffs_keep));

for i=1:size(c2r,1)
    num_null = 300 - sum(isnan(c2r_null_coeffs_keep(i,:)));
    c2r(i,:) = sum(c2r_null_coeffs_keep(i,:) > best_c2r_coeffs_keep(i)) ./ num_null;
    num_null = 300 - sum(isnan(c2r_null_coeffs_keep(i,:)));
    r2c(i,:) = sum(c2r_null_coeffs_keep(i,:) > best_r2c_coeffs_keep(i)) ./ num_null;
end

[~, prod_sort] = sort(mean_FR_keep(:,1) .* mean_FR_keep(:,2));
[~, RFA_sort ] = sort(mean_FR_keep(:,2));
[~, CFA_sort ] = sort(mean_FR_keep(:,1));
[~, low_sort] = sort(min(mean_FR_keep, [],2));
[~, which_region] = min(mean_FR_keep(low_sort, :), [], 2);

%%

%%
close all

fig1 = figure;
%{
total_pairs = size(low_sort);
num_keeps = floor(low_sort/2);
num_per_region = floor(num_keeps/2);

temp1 = low_sort(which_region==1);
temp2 = low_sort(which_region==2);

keep = [temp1(num_per_region:end); temp2(num_per_region:end)];
%}

%last = size(this_sort, 1);

last = size(low_sort, 1);
start = floor(.5 * last);
keep = low_sort(start:last);

final_c2r_coeffs = best_c2r_coeffs_keep(keep);
final_r2c_coeffs = best_r2c_coeffs_keep(keep);

c2r = c2r(keep);
r2c = r2c(keep);

histogram(c2r, 'Facecolor', 'blue', 'Facealpha', .5, 'BinWidth', 0.02, 'Normalization','probability')
xlim([0,1])
hold on
histogram(r2c, 'Facecolor', 'red', 'Facealpha', .5, 'BinWidth', 0.02, 'Normalization','probability')
xlim([0,1])
title('INTER P-VALUES, -300->0: blue CFA->RFA, red RFA->CFA')
saveas(fig1, 'C:\Users\mirilab\OneDrive - Northwestern University\Desktop\diya\CCM\exclusion testing\low sort\p.png')

%%
c2r_coeffs = final_c2r_coeffs;
r2c_coeffs = final_r2c_coeffs;
%%
fig2 = figure
histogram(c2r_coeffs, 'Facecolor', 'blue', 'Facealpha', .5, 'BinWidth', 0.02)
xlim([0,1])
hold on
histogram(r2c_coeffs, 'Facecolor', 'red', 'Facealpha', .5, 'BinWidth', 0.02)
xlim([0,1])
title('INTER STRENGTHS: blue CFA->RFA, red RFA->CFA, exclude low of pair FRs')
saveas(fig2, 'C:\Users\mirilab\OneDrive - Northwestern University\Desktop\diya\CCM\exclusion testing\low sort\str.png')


%% FRS
fig3= figure
CFA_FR_exclude = mean_FR_keep(keep, 1);
RFA_FR_exclude = mean_FR_keep(keep, 2);
histogram(log10(CFA_FR_exclude), 'Facecolor', 'blue', 'Facealpha', .5, 'BinWidth', 0.05)
hold on
histogram(log10(RFA_FR_exclude), 'Facecolor', 'red', 'Facealpha', .5, 'BinWidth', 0.05)
title('CFA=Blue, RFA=red, exclude low of pair FRs')
xlabel('log10 FR')
saveas(fig3, 'C:\Users\mirilab\OneDrive - Northwestern University\Desktop\diya\CCM\exclusion testing\low sort\fr.png')



