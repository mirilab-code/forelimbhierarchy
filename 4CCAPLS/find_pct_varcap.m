function [meanpctVC, sempctVC] = find_pct_varcap(pcs,analyzed)


cs_varcap_pcs = cumsum(pcs,2);
cs_varcap_test = cumsum(analyzed,2);

% get the 10th PC
pc10_vc_pcs = cs_varcap_pcs(:,10);
pc10_vc_test = cs_varcap_test(:,10);

pct_varcap = pc10_vc_test ./ pc10_vc_pcs;
meanpctVC = mean(pct_varcap);
stdpctVC = std(pct_varcap);
sempctVC = stdpctVC / sqrt(length(pct_varcap));

fprintf('%0.3f +/- %0.3f \n',meanpctVC,sempctVC);

end