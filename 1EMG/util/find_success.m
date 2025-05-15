function success_bounds = find_success(num_spouts,reward_times,touch,suction,err)
    % INPUT: the channel#, the reward times, all the touch
    % channels, and the suction channel. If the channel of the first touch 
    % AFTER the reward channel is the same as the channel# i.e. if spout 3
    % is touched and the channel# is 3, AND this happens before the suction
    % channel, then it's a success.
    %
    % OUTPUT: the #successes x 2 matrix: each row being [time of reward, time of touch]
    
    success_bounds = cell(1,num_spouts);
    for ii = 1:num_spouts
        channel = ii;
        touch = [touch; suction];
        nChannels = size(touch,1);

        for jj=1:length(reward_times{ii})
            t = reward_times{ii}(jj);
            suc = suction(t:end);
            stop = find(suc,1);
            if (isempty(stop))
                stop = length(touch(channel,t:end))-1;
            end

    %         disp([t t+stop]);
            trace = touch(:,t:t+stop);

            reaches = zeros(nChannels,1);
            for c=1:nChannels
                q = find(trace(c,:)~=0,1);
                if (~isempty(q))
                    reaches(c) = q;      
                else
                    reaches(c) = NaN;
                end
            end

            % Check if the first reach is at the correct channel
            % Also check if the correct spout is touched within a few ms 
            %  of an incorrect spout
            first_touch = min(reaches);
            first_reach = find(reaches==first_touch);

    %         disp(reaches);
            correct_touch = reaches(channel);
            diff_touches = correct_touch-first_touch;
    %         disp(diff_touches);
            if (diff_touches <= err)
    %             disp('success!');
                bounds = [t t+correct_touch];
                success_bounds{ii} = [success_bounds{ii}; bounds];
            end
        end
    
    end

end
