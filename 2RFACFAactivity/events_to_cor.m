function [CFA,RFA,u1,u2]= events_to_cor(events, t1, outlier_thresh)

%INPUTS events and t for threshold (for now, 1500)

% LINES 34 and 37 NEED MANUAL INPUTS 


%OUTPUTS CFA and RFA trains that can be run through JCCG
%need to check that u1 and u2 match with histogram

e1 = events{1}; %pulls out CFA
u1= max(e1(:,3))- t1; %sets limit for CFA 

cortical_neurons1 = find(e1(:,3)>=u1);% finds cortical neurons
st1 = e1(cortical_neurons1,2);
index1= e1(cortical_neurons1,1);


A= [index1'; st1'];
CFA_ecor=A';

CFA= events_to_train(CFA_ecor); % gets spike trains for CFA




%gets rid of outliers
bad_inds= find(events{2}(:,3)>outlier_thresh); 
events{2}(bad_inds, :)= [];

e2 = events{2}; % pulls out RFA

edges =[0:1:max(e2(:,3))];% if no manual input, make it max of e2 (:,3), otherwise same as line 34
H= histcounts (e2(:,3),edges); %makes hist of spike depth

% finding region of 200 or more 0s in spike depth histo
props = regionprops(logical(H ==0), 'Area', 'PixelIdxList');
indexesOf200orMore = find([props.Area]>= 200);
for k = indexesOf200orMore
  theseIndexes = props(k).PixelIdxList;
  A= (theseIndexes); %array of indexes with 200 0s in a row
end

u2= max(A)-100; %pulls out max index of gap

if length(u2)>1
    u2=0;
end

if u2==0
    u2= max(e2(:,3))-t1;
end

cortical_neurons2 = find(e2(:,3)>=u2);
st2 = e2(cortical_neurons2,2);
index2= e2(cortical_neurons2,1);


B= [index2'; st2'];
RFA_ecor=B';
RFA= events_to_train(RFA_ecor);

%u2 %reports what threshold was used for RFA 

end


