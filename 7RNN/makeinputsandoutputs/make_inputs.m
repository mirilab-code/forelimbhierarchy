%% create the go cue and laser signals
clc


% make the laser signal. starts at 0, ramps up for 15ms so that it reaches
% 1 by 15ms after onset. decays quicker after 50ms.

sigmoid = @(x,slope,shift) 1./(1 + exp(-slope.*(x-shift)));
x = -50:50;
beginning = sigmoid(x,.5,8);


startramp = find(beginning>0.01,1)
topout = find(beginning>0.95,1)
w = topout - startramp;

endx = 1:50;
endpart = 1-sigmoid(endx,1,0);

y = [beginning endpart];

figure;
hold on
plot(y)
xline(50)
xline(startramp, '--r')
xline(topout, '--r')
xline(100)
xline(105)
hold off

% when inactivating RFA just send 0s to CFA and vice versa.
input.CFAinact.cfa = y;
input.CFAinact.rfa = y*0;
input.RFAinact.cfa = y*0;
input.RFAinact.rfa = y;


%% make the go cue
clc

%% 25ms before laser
x = -50:100;
go25 = sigmoid(x,0.3,-26);

startramp = find(go25>0.01,1)
half = 50-find(go25==0.5,1);
fprintf('%dms before the laser \n',half)
topout = find(go25>0.95,1)

figure;
hold on
plot(go25)
xline(50)
xline(startramp,'--r')
xline(50-half, '--r')
xline(topout, '--r')
hold off 

%% 20ms before laser
x = -50:100;
go20 = sigmoid(x,0.3,-21);

half = 50-find(go20==0.5,1);
fprintf('%dms before the laser \n',half)
% topout = find(y>0.95,1)

figure;
hold on
plot(go20)
xline(50)
xline(50-half, '--r')
% xline(topout, '--r')
hold off

%% 15ms before laser
x = -50:100;
go15 = sigmoid(x,0.3,-16);

half = 50-find(go15==0.5,1);
fprintf('%dms before the laser \n',half)
% topout = find(y>0.95,1)

figure;
hold on
plot(go15)
xline(50)
xline(50-half, '--r')
% xline(topout, '--r')
hold off

%% 10ms before laser
x = -50:100;
go10 = sigmoid(x,0.3,-11);

half = 50-find(go10==0.5,1);
fprintf('%dms before the laser \n',half)
% topout = find(y>0.95,1)

figure;
hold on
plot(go10)
xline(50)
xline(50-half, '--r')
% xline(topout, '--r')
hold off



%%
input.go25 = go25;
input.go20 = go20;
input.go15 = go15;
input.go10 = go10;

save('model_setup_data\input.mat','input');
disp('saved!')












%% 
disp('done!')