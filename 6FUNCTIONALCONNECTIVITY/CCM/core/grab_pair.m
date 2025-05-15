function [CFA_FR, RFA_FR, reaches] = grab_pair(row_top10k)
%row top10k is a single row from top_ten_thou_upd structure

folder_list = load('./10k_data/folder_names.mat').output;

folder_number = row_top10k(1);
CFA_unit = row_top10k(3);
RFA_unit = row_top10k(4);

date = folder_list(folder_number, :);
temp = load(sprintf('./10k_data/%s/%s_pyr_CFA_RFA_trains_sort.mat', date, date)).pyr_CFA_RFA_trains_sort;

duration = max(cell2mat(temp{1,1}));

reaches = load(sprintf('./10k_data/%s/%s_reach_bounds_edit.mat', date, date)).reach_bounds_edit;

CFA_FR = singletrain_to_firingrate(temp{1,1}{CFA_unit, 1}, duration);
RFA_FR = singletrain_to_firingrate(temp{1,2}{RFA_unit, 1}, duration);

end

