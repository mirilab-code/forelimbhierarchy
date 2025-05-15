function [x_FR, y_FR, reaches] = grab_pair_intra(row_intrapairs)
%row top10k is a single row from top_ten_thou_upd structure

folder_list = load('./10k_data/folder_names.mat').output;

folder_number = row_intrapairs(1);
x = row_intrapairs(2);
y = row_intrapairs(3);
region = row_intrapairs(4);

date = folder_list(folder_number, :);
temp = load(sprintf('./10k_data/%s/%s_pyr_CFA_RFA_trains_sort.mat', date, date)).pyr_CFA_RFA_trains_sort;

duration = max(cell2mat(temp{1,1}));

reaches = load(sprintf('./10k_data/%s/%s_reach_bounds_edit.mat', date, date)).reach_bounds_edit;

x_FR = singletrain_to_firingrate(temp{1,region}{x, 1}, duration);
y_FR = singletrain_to_firingrate(temp{1,region}{y, 1}, duration);

end

