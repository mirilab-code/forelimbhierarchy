function waveform_dir = load_waveformDirection(date,animal,units)
%load_eaveformDirection calls the waveformDirection function.
%   date: is an 8 character array of the date of recording
%   animal: is a character array of the animal name
%   units: is a cell array of all the units being analyzed (must be a cell
%   array of length 2)

addpath('Z:\Scripts\neuropixel');
kilosort_dir = 'kilosort';
waveform_dir  = cell(1,length(units));
for ii = 1:length(units)
str = sprintf('Z:\\neuropixel\\newMaps\\%simec%d_neuropixPhase3B1_kilosortChanMap.mat',...
    animal,ii-1);
channelMapFile = str;
imec_path = sprintf('Z:\\akiko\\%s_%s_g0\\%s_%s_g0_imec%d\\%s_%s_g0_t0.imec.ap.bin',...
    date,animal,date,animal,ii-1,date,animal);
ks_path = sprintf('Z:\\akiko\\%s_%s_g0\\%s_%s_g0_imec%d\\%s',...
    date,animal,date,animal,ii-1,kilosort_dir);
waveform_dir{ii} = waveformDirection(channelMapFile,...
    imec_path,ks_path,units{ii});
end

