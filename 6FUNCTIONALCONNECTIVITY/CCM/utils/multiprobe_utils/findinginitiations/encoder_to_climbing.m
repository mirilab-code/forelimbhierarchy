function climbing_bouts = encoder_to_climbing(analogin, maxT)
%
%
%
%
addpath(genpath('Z:/code'));
%% diya modifications of natalies code. i use this for TE
%get encoder signal
enc_ch = 2;

encoderSignal = analogin(enc_ch, 1:maxT)'; 

% Wheel velocity (encoder) per session
encoder = rescale(encoderSignal); % rescale 0 to 1 (no longer in pure V)
%vel_boxcar_span_samps = vel_boxcar_span * sampleRate_ds / 1000;
signal = smoothdata(encoder, 1, 'movmedian', 100); %not sure if i need this
signal = smoothdata(signal, 1, 'movmean', 5000); 
M = movstd(signal, 10000);

%figure; plot(M); hold on; yline(0.005)

%%
thr = 0.005;
M2 = find(M >= thr); %these are all the moments where moving
M3_end = find(diff(M2) > 1); % end of bouts, whenever theres a large amount of time without movement
M3_init = M3_end+1;
M3idx_end = [M2(M3_end); M2(end)]; 
M3idx_init = [M2(1); M2(M3_init)]; 
if length(M3idx_init)>length(M3idx_end); M3idx_init(end) = [];end
stitch_bouts = [M3idx_init, M3idx_end, M3idx_end-M3idx_init];

% Get the difference between each start/end pt
thresh_diff = 0.01;
for k = 1:length(stitch_bouts)
   stitch_bouts(k,4) = mean(M(stitch_bouts(k,1):stitch_bouts(k,2))); 
end
t = find(stitch_bouts(:,4) < thresh_diff);
stitch_bouts(t, :) = [];

%figure; plot(signal); hold on;
%plot(stitch_bouts(:,1), signal(stitch_bouts(:,1)),  'b*')
%plot(stitch_bouts(:,2), signal(stitch_bouts(:,2)),  'r*')
%title([animal ' sess' num2str(i)]); 
%legend('signal', 'start bout', 'end bout')

%% 
% if the signal change is negative, we need to split a bout
% so that we can calculate the accurate dist covered later
for k = 1:length(stitch_bouts)
    stitch_bouts(k,7) = signal(stitch_bouts(k,1)) - signal(stitch_bouts(k,2));
    if stitch_bouts(k,7) < 0
       smin = min(signal( stitch_bouts(k,1):stitch_bouts(k,2)));
       smax = max(signal( stitch_bouts(k,1):stitch_bouts(k,2)));
       sminidx = find(signal( stitch_bouts(k,1):stitch_bouts(k,2)) == smin) + stitch_bouts(k,1) - 1;
       smaxidx = find(signal( stitch_bouts(k,1):stitch_bouts(k,2)) == smax) + stitch_bouts(k,1) - 1;
       stitch_bouts(k,8) = min(sminidx);
       stitch_bouts(k,9) = min(smaxidx);
    end
end

% stitch_bouts matrix: 
% 1. start bout idx
% 2. end bout idx
% 3. duration of bout in samples
% 4. difference in movstd between start and end of bout
% 5. distance covered in bout (cm)
% 6. duration of bout in seconds
% 7. difference in encoder signal between start and end of bout
% 8. idx of min point of bout (only when encoder hits floor and starts
% at the top)
% 9. idx of max point of bout

%% probably dirty and inefficient but lets convert stitch_bouts to logical array
climbing_bouts = zeros(length(encoderSignal), 1, 'logical');
for i=1:length(stitch_bouts)
    climbing_bouts(stitch_bouts(i, 1):stitch_bouts(i,2)) = logical(1);
end