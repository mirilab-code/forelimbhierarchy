mean_FR = zeros(6020,2);
pb = progressbar(' getting firing rate... ');
for i=1:size(intra_pairs, 1)
    pb.print(i, size(intra_pairs, 1));
    [X_FR, Y_FR, reaches] = grab_pair_intra(intra_pairs(i,:));
    [~, bouts] = get_reaching_bouts(reaches, size(X_FR,2), -200, 799);
    mean_FR(i, 1) = mean(X_FR(bouts));
    mean_FR(i, 2) = mean(Y_FR(bouts));
end
save('mean_FR_intra.mat', 'mean_FR')

temp_CFA_coeffs = [best_x2y_coeffs(1:3010); best_y2x_coeffs(1:3010)];
temp_RFA_coeffs = [best_x2y_coeffs(3011:end); best_y2x_coeffs(3011:end)];

temp_CFA_nulls = [x2y_null_coeffs(1:3010, :); y2x_null_coeffs(1:3010, :)];
temp_RFA_nulls = [x2y_null_coeffs(3011:end, :); y2x_null_coeffs(3011:end, :)];

CFA_mean_FR = [mean_FR(1:3010, 1); mean_FR(1:3010, 2)];
RFA_mean_FR = [mean_FR(3011:end, 1); mean_FR(3011:end, 2)];

best_c2c_coeffs_keep = temp_CFA_coeffs(~isnan(temp_CFA_coeffs));
CFA_mean_FR_keep = CFA_mean_FR(~isnan(temp_CFA_coeffs));
CFA_nulls_keep = temp_CFA_nulls(~isnan(temp_CFA_coeffs), :);

best_r2r_coeffs_keep = temp_RFA_coeffs(~isnan(temp_RFA_coeffs));
RFA_mean_FR_keep = RFA_mean_FR(~isnan(temp_RFA_coeffs));
RFA_nulls_keep = temp_RFA_nulls(~isnan(temp_RFA_coeffs), :);

%[~, prod_sort] = sort(mean_FR_keep(:,1) .* mean_FR_keep(:,2));
[~, RFA_sort ] = sort(RFA_mean_FR_keep);
[~, CFA_sort ] = sort(CFA_mean_FR_keep);
%[~, low_sort] = sort(min(mean_FR_keep, [],2));


c2c = zeros(size(best_c2c_coeffs_keep));
r2r = zeros(size(best_r2r_coeffs_keep));
for i=1:size(c2c,1)
    num_null = 300 - sum(isnan(CFA_nulls_keep(i,:)));
    c2c(i,:) = sum(CFA_nulls_keep(i,:) > best_c2c_coeffs_keep(i)) ./ num_null;
end

for i=1:size(r2r, 1)
    num_null = 300 - sum(isnan(RFA_nulls_keep(i,:)));
    r2r(i,:) = sum(RFA_nulls_keep(i,:) > best_r2r_coeffs_keep(i)) ./ num_null;
end

CFA_end = size(CFA_sort,1);
RFA_end = size(RFA_sort,1);

CFA_start = floor(CFA_end/2);
RFA_start = floor(RFA_end/2);

c2c_coeffs = best_c2c_coeffs_keep(CFA_sort(CFA_start:CFA_end));
r2r_coeffs = best_r2r_coeffs_keep(RFA_sort(RFA_start:RFA_end));

c2c = c2c(CFA_sort(CFA_start:CFA_end));
r2r = r2r(RFA_sort(RFA_start:RFA_end));

%c2r = best_c2c_coeffs_keep(this_sort(start:last));
%r2c = best_r2r_coeffs_keep(this_sort(start:last));

histogram(c2c, 'Facecolor', 'blue', 'Facealpha', .5, 'BinWidth', 0.05)
xlim([0,1])
hold on
histogram(r2r, 'Facecolor', 'red', 'Facealpha', .5, 'BinWidth', 0.05)
xlim([0,1])
title('INTRA: blue CFA->CFA, red RFA->RFA, sort by lower FR of intra pair, exclude bottom 50%')


