%% analyze the new and old weights

oldW = readNPY('oldW.npy');
newW = readNPY('newW.npy');


%%
diffW = newW - oldW;

imagesc(diffW);





