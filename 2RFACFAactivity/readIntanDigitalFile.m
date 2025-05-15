function dig_in = readIntanDigitalFile(fname, num_channels)
%
% function to read in a digital IO Intan file ('digitalIn' or 'digitalOut')
%
% INPUTS
%   fname - filename of the Intan digital file to read
%
% OUTPUTS
%   digital_word - vector of uint16's containing bit-wise values of each
%      digital line from the Intan system

fileinfo = dir(fname);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(fname, 'r');
dig_in_raw = fread(fid, num_samples, 'uint16');
dig_in = zeros(num_channels, num_samples);
fclose(fid);

for i=1:num_channels
   mask = 2^(i) * ones(size(dig_in_raw));
   dig_in(i, :) = (bitand(dig_in_raw, mask) > 0);
end


end