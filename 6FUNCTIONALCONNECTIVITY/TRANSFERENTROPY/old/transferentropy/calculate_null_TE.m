function [obsTE, nullTE, shifts] = calculate_null_TE(stacked_bintrain,times,delay,iter, calcObs)

% stacked_bintrain: stacked binary train where each train is the
%                   concatenated epochs we want to analyze.
% times: a vector contained where to insert the padding in each binary train
% delay: how much to delay to calcuate for the TE which is also how much
%        padding of 0z to add between epochs
% iter: # of null maps to make
% calcObs: whether or not to calculate the observed TE

tic

N = size(stacked_bintrain,1);


% get the observed TE first
if(calcObs)
    disp('calculating the observed transfer entropy...');

    paddedBinTrain = [];
    for i=1:N
        pbt = insert_padding(stacked_bintrain(i,:),times,delay);
        paddedBinTrain = [paddedBinTrain; pbt];
    end

    [obsTE,~,~] = calculateTE(paddedBinTrain,delay);    
    
else 
    obsTE = 0;
end

% now do all the null map stuff
nullTE = zeros(N,N,iter);
shifts = logical(nullTE*0);

disp('creating null map...');
for i=1:iter
    disp([i iter]);
    [shiftedTE,allowed_shifts] = calculate_shifted_TE(stacked_bintrain,times,delay);
    nullTE(:,:,i) = shiftedTE;
    shifts(:,:,i) = allowed_shifts;
end

toc

end