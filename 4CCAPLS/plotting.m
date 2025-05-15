%% plotting
hold on
for i=1:size(CFA_RG,3)
    A = CFA_RG(:,:,i);
    A = mean(A,1);
    if(i==61)
        plot(A, 'LineWidth',1, 'Color', [0, 0, 1, 1])
    else
        plot(A, 'LineWidth',1, 'Color', [0.4, 0.4, 0.4, 0.5])
    end
    
end

yyaxis right
B = mean(RFA_RG,1);
plot(B, 'LineWidth',1,'Color',[1 0 0 1])


%%