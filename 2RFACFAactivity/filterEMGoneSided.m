function [ out ] = filterEMGoneSided( in , varargin)
% An updated version of filterEMGoneSided/filterEMG. Now change
% samplingRate also as an input variable since the rate of in is not always
% 1000 Hz. 
% filters EMG with three steps: 1) high pass; 2) rectify; 3) low pass 

% Please change: filterEMG(EMG_traces(i,:),high_pass_f,std_gauss_scaled),
% or filterEMGoneSided(EMG_traces(i,:),high_pass_f,std_gauss_scaled) to
% filterEMGoneSided_updated(EMG_traces(i,:),high_pass_f,std_gauss_scaled,sample_rate)
%% Set initial parameters and parse command line
if nargin ==4;
    hp = varargin{1};
    SD = varargin{2};
    samplingRate = varargin{3};
elseif nargin == 3;
    hp = varargin{1};
    SD = varargin{2};
    samplingRate = 1000;
elseif nargin == 2;
    hp = varargin{1};
    SD = 25; % default low pass to 25;
    samplingRate = 1000;
elseif nargin == 1;
    hp = 40;
    SD = 25;
    samplingRate = 1000;
else
    error('Unrecognized number of arguments')
end

filterOrder = 12; % we have tried a bunch of different filter orders, doesn't make much difference

%% 1) High pass: design a high pass butterworth filter with cutoff frequency of hp
d = fdesign.highpass('N,F3db',filterOrder,hp,samplingRate);
emg_hp = design(d, 'butter');
in_hp = filtfilt(emg_hp.sosMatrix, emg_hp.ScaleValues, in);

%% 2) Rectify
in_rec = abs(in_hp); %recitfy

%% 3) Low Pass
out = filterGaussOneSided(in_rec, SD);

end

