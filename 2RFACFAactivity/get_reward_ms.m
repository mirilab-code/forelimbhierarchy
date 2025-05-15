



% filebase = 'a041_191218_111855';



function get_reward_ms (dir_base)

% if nargin == 1
%     dir_base = 'Z:\akiko/EMG/';
%     
% end
% 
% cd (dir_base);
% 
%%

overwrite = 1;
S = get(0, 'ScreenSize');

%%
% if ~exist('fig','dir')
%     mkdir('fig');
% end
% cd ([dir_base,'/fig']);
% cd ../
% 
cd (dir_base);
if ~exist('event','dir')
    mkdir('event');
end
if ~exist('trace_each','dir')
    mkdir('trace_each');
end


%%
cd (dir_base);
if overwrite == 1 
    % import     
    [frequency_parameters,board_adc_channels,~]=read_Intan_RHD2000_file_extract_auto;% read info.rhd file from current folder automatically
    % sample rate
    r = frequency_parameters.board_adc_sample_rate; % sample rate from header file
    r = r./1000;%kHz
    fileinfo = dir('analogin.dat');
    if isempty(fileinfo)
        fprintf('ERROR; There are no analogin.dat file.\n');
        return;
    end
    %num_channels = length(board_adc_channels); % ADC input info from header file, 1)wheel, 2)RW, 3)Tone, 4)Brake, 5)LED
    
    if ~exist('digitalin.dat','file')
        fprintf('digitalin.dat does not exist!\n');
    else
        num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
        fid = fopen('digitalin.dat', 'r');
        digital_word = fread(fid, num_samples, 'uint16');% digital outputs channel is 
        fclose(fid);
        %
        %  value. 
        %{  
        Ch0: none
        Ch1: Touch sensor 4
        Ch2: Touch sensor 3
        Ch3: Touch sensor 2
        Ch4: Touch sensor 1
        Ch5: Resting bar touch sensor
        Ch6: Reward 1 => 2^6 = 64
        Ch7: Reward 2 => 2^7 = 128
        Ch8: Reward 3 => 2^8 = 256
        Ch9: Reward 4 => 2^9 = 512
        Ch10: Suction
        Ch11-15: none        
        %}
        
        % no down sample,
        % sampling rate is still 20kHz
        alltime = [1:length(digital_word)]';
        d_rw1only = alltime(digital_word == 64,1);% two value around the rw1, 64 and 96, 96-64=32 = 2^5, Ch5 touch sensor(only several ms).
        d_rw2only = alltime(digital_word == 128,1);
        d_rw3only = alltime(digital_word == 256,1);
        d_rw4only = alltime(digital_word == 512,1);
        
        % find the event 
        % reward signal must continues ~200ms = 200x20 = 4000 .. <= 2000 to
        % 6000, around 4000/20 = 200ms. OK
        
        % delete the event within 500ms (10000pts)
        del_pts = 10000;
        
        del_event = [];
        for i = 2:size(d_rw1only,1)
            if d_rw1only(i) - d_rw1only(i-1) < del_pts
                del_event = [del_event;i];
            end            
        end
        d_rw1_start = d_rw1only;
        d_rw1_start(del_event,:) = [];% still dat pts
        if isempty(d_rw1_start)
            fprintf('reward1 does not exist!\n');
        end
        
        del_event = [];
        for i = 2:size(d_rw2only,1)
            if d_rw2only(i) - d_rw2only(i-1) < del_pts
                del_event = [del_event;i];
            end            
        end
        d_rw2_start = d_rw2only;
        d_rw2_start(del_event,:) = [];% still dat pts
        if isempty(d_rw2_start)
            fprintf('reward2 does not exist!\n');
        end
        
        del_event = [];
        for i = 2:size(d_rw3only,1)
            if d_rw3only(i) - d_rw3only(i-1) < del_pts
                del_event = [del_event;i];
            end            
        end
        d_rw3_start = d_rw3only;
        d_rw3_start(del_event,:) = [];% still dat pts
        if isempty(d_rw3_start)
            fprintf('reward3 does not exist!\n');
        end
               
        del_event = [];
        for i = 2:size(d_rw4only,1)
            if d_rw4only(i) - d_rw4only(i-1) < del_pts
                del_event = [del_event;i];
            end            
        end
        d_rw4_start = d_rw4only;
        d_rw4_start(del_event,:) = [];% still dat pts
        if isempty(d_rw4_start)
            fprintf('reward4 does not exist!\n');
        end
        clear('d_rw1only','d_rw2only','d_rw3only','d_rw4only');
        
        total_rw = length(d_rw1_start) + length(d_rw2_start) + length(d_rw3_start) + length(d_rw4_start)
        
        % dat pts to ms
        rw1 = round(d_rw1_start./r);
        rw2 = round(d_rw2_start./r);
        rw3 = round(d_rw3_start./r);
        rw4 = round(d_rw4_start./r);
        
        cd ([dir_base '/event']);
        save('rw_ms', 'rw*');% sampling rate is 1kHz
       
    end
    
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    