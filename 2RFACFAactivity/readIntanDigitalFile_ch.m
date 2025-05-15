function dig_in = readIntanDigitalFile_ch(fname, board_dig_in_channels)
%
% function to read in a digital IO Intan file ('digitalIn' or 'digitalOut')
%
% INPUTS
%   fname - filename of the Intan digital file to read
%   board_dig_in_channels - from read_Intan_RHD2000_file_extract_auto2
%   

% OUTPUTS
%   digital_word - vector of uint16's containing bit-wise values of each
%      digital line from the Intan system

num_channels = length(board_dig_in_channels);
ch_num = zeros(1,num_channels);
for i=1:num_channels
     ch_num(i) = board_dig_in_channels(i).native_order;
end

fileinfo = dir(fname);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(fname, 'r');
dig_in_raw = fread(fid, num_samples, 'uint16');
dig_in = zeros(num_channels, num_samples);
fclose(fid);

row_num = 0;
for i=ch_num
   row_num = row_num+1; 
   mask = 2^(i) * ones(size(dig_in_raw));
   dig_in(row_num, :) = (bitand(dig_in_raw, mask) > 0);    
end




end