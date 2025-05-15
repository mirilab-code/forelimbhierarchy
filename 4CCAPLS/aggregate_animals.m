function M = aggregate_animals(D)

Da{1} = D(1:6);
Da{2} = D(7:8);
Da{3} = D(9:13);
Da{4} = D(14:17);
Da{5} = D(18:19);
Da{6} = D(20:end);

%%
for i=1:length(Da)
    disp(i);
    d = Da{i};
    
    tr = cellfun(@(x) size(x,2), {d.RFA_RG});
    tc = cellfun(@(x) size(x,2), {d.CFA_RG});
    t = [tc;tr]';
    exc = mean(t,2)~=1600;
    d(exc) = [];


    x = {d.RFA_RG};
    RFA_RG = cat(1,x{:});
    x = {d.CFA_RG};
    CFA_RG = cat(1,x{:});

    m.RFA_RG = RFA_RG;
    m.CFA_RG = CFA_RG;
    m.lags = d.lags;

    if(i==1)
        M = m;
    else
        M = [M m];
    end
end

%%
% D = M;


end