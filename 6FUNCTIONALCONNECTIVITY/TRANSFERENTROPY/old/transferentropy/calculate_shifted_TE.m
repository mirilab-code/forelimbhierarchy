function [TEshift,effective_shifts] = calculate_shifted_TE(stacked_bintrain,times,delay)

% tic

N = size(stacked_bintrain,1);
T = size(stacked_bintrain,2);
shifts = 3000:T-3000;   % acceptable shifts are more than 3seconds or less than T-3 seconds
unit_shifts = randsample(shifts,N,true);    % true = with replacement

effective_shifts = zeros(N,N);

catBinary = [];
for i=1:N
    shift = unit_shifts(i);
    shiftedtrain = circshift(stacked_bintrain(i,:),shift);
    catshift = insert_padding(shiftedtrain,times,delay);
    catBinary = [catBinary; catshift];
end

for i=1:N
    for j=1:N
        effective_shifts(i,j) = abs(unit_shifts(i)-unit_shifts(j));
    end
end

effective_shifts = (effective_shifts>3000 & effective_shifts<(T-3000));

% now calculate the TE with the shifted trains
[TEshift,~,~] = calculateTE(catBinary,delay);

% toc
end