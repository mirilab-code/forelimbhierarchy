
D_vc = {};
D_vc{1} = diff(a048_vc_CFARFA');
D_vc{2} = diff(a050_vc_CFARFA');
D_vc{3} = diff(a051_vc_CFARFA');
D_vc{4} = diff(MA1_vc_CFARFA');
D_vc{5} = diff(MA2_vc_CFARFA');
D_vc{6} = diff(ss2_vc_CFARFA');

hold on
% cellfun(@plot, D_vc);
cellfun(@(B) plot(B, 'ok'), D_vc)
xlim([0 10])

disp('done!')