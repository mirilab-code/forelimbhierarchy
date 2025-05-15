
P_values= [];
allowed_shifts= [];
TE_nulls = [];
CtoCs= [];


for i= 1: length(sess)
    
    disp('%i')
    P=[];
    all_shifts=[];
    TE_null= [];
    CtoC=[];
    
    [CtoC,TE_null, all_shifts]= TE_analysis_HS (sess{i},CFApairs, i);
     P = get_pvals(CtoC,TE_null,all_shifts);
    
    P_values= [P_values; P];
    allowed_shifts= [allowed_shifts; all_shifts];
    TE_nulls= [TE_nulls; TE_null];
    CtoCs= [CtoCs; CtoC];
    
end