function movement= movement_only(EMG,long_window_cutoff,short_window_cutoff)

EMG_sum= sum(EMG);



%have to identify stretch of inactivity first (3005000 - 3010000 for 0327)
lb = input('not moving segment start sample: '); 
ub = input('not moving segment end sample: ');
sd = std(EMG_sum(lb:ub));
avg= mean(EMG_sum(lb:ub));
thresh = avg+(7*sd); 


%classfication
movement = EMG_sum > thresh; 

%readding in patches of zeros that are less than long_window_cutoff ms
props = regionprops(logical(movement==0), 'Area', 'PixelIdxList');
indexesOf100orless = find([props.Area]<= long_window_cutoff);
for k = indexesOf100orless
  theseIndexes = props(k).PixelIdxList;
  movement(theseIndexes) = 1;
end
 
%getting rid of patches too small to be movement
props = regionprops(logical(movement==1), 'Area', 'PixelIdxList');
indexesOf10orless = find([props.Area]<= short_window_cutoff);
for k = indexesOf10orless
  theseIndexes = props(k).PixelIdxList;
  movement(theseIndexes) = 0;
end

movement = logical(movement);

%sanity check
plot(EMG_sum);
hold on
plot(movement*100);
end



