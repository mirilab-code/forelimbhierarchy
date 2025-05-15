%% make model data
load('C:\Users\mirilab\OneDrive - Northwestern University\Documents\Mark\RNN_rfacfa\model_setup_data\input.mat')
load('C:\Users\mirilab\OneDrive - Northwestern University\Documents\Mark\RNN_rfacfa\model_setup_data\network.mat')
load('C:\Users\mirilab\OneDrive - Northwestern University\Documents\Mark\RNN_rfacfa\model_setup_data\output.mat')
window = -50:100;
save('model_setup_data\data.mat','input','output','network','window');
disp('done making data.mat')