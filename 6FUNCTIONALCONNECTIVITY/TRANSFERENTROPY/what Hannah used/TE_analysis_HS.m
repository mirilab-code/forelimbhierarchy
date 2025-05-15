
function [RtoC,TE_null, all_shifts]= TE_analysis_HS (binary_trains_all_trials, top_ten_thou_upd, rec_num)

%%inputs- binary_trains_all_trails and top_ten_tho_upd from shared drive,
%%rec_num is the recording number based on top_ten_tho_upd  

%%outputs- RtoC is number of pairs x 3 matrix. first two columns are RFA
%%and CFA neuron index respectively. third column is TE value
%%TE_null is number of pairs x 1 matrix. Each elements has all 300 values
%%of nullTE
%%allowed_shifts is


%I SWAPPED RFA AND CFA, PLEASE CHECK IT WHEN IT IS DONE RUNNING
tic

dimm = size(binary_trains_all_trials{1});
nRFA= dimm(1);
dim = size(binary_trains_all_trials{1});
nCFA= dim(1);
nUnits = nRFA + nCFA; 


stacked_initiation= [];

for i= 1: nRFA
    S =binary_trains_all_trials{2} (i, :,:);
    S=reshape(S, [1, dimm(2)*dimm(3)]);
    stacked_initiations (i,:)=S;
end 
    
for i= 1: nCFA
    S =binary_trains_all_trials{1} (i, :,:);
    S=reshape(S, [1, dim(2)*dim(3)]);
    stacked_initiations (i+nRFA,:)=S;
end 

  slices= [1:1001:length(stacked_initiations)];

F = find(top_ten_thou_upd(:,1)==rec_num);
pairs= [];

for i= 1: length(F) %SWAPPED FOR RFA?CFA SWItch
pairs(i,1) = top_ten_thou_upd(F(i),2);
pairs(i,2) = top_ten_thou_upd(F(i),3); 
end
  
toc
%%  testing
tic
[observedTE,nullTE,allowed_shifts] = calculate_null_TE(stacked_initiations,slices,30,15,true);


%  testing 2
[~,nullTE2,allowed_shifts2] = calculate_null_TE(stacked_initiations,slices,30,285,false);

nullTE = cat(3,nullTE,nullTE2);
allowed_shifts = cat(3,allowed_shifts,allowed_shifts2);

disp ('now we find pairs')



for i = 1:length(pairs)
    TE_final (i) = observedTE (pairs(i,1) , pairs(i,2)+nRFA );
     b = nullTE (pairs(i,1) , pairs(i,2)+nRFA, : );
     c= allowed_shifts (pairs(i,1),pairs(i,2)+nRFA, :);
     TE_null{i,1}= reshape(b, [1, length(b)]);
     all_shifts{i,1} = reshape(c,[1,length(c)]);
end

RtoC= [pairs, TE_final'];

toc
end




%%