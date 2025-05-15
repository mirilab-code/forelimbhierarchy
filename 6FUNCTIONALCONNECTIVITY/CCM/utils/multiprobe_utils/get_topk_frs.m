function [RFA_idx, CFA_idx] = get_topk_frs(RFA_FRs,CFA_FRs, k)
%GET_TOP10K_FR Summary of this function goes here
%   Detailed explanation goes here
%k=50; %set how many top

samples = size(RFA_FRs, 2);
fr_pairs = mean(RFA_FRs, 2) * mean(CFA_FRs, 2)';
[~, lin_idx] = maxk(fr_pairs(:), k);
[RFA_idx, CFA_idx] = ind2sub(size(fr_pairs), lin_idx);


end
