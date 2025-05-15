function events = removeFlaggedUnits(events,unit_flags)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
num_probes = length(events);
events_flags = cell(num_probes,1);
for ii = 1:num_probes
    for jj = 1:length(events{ii}(:,1))
        for kk = 1:length(unit_flags{ii})
            if events{ii}(jj,1) == unit_flags{ii}(kk)
                events_flags{ii} = [events_flags{ii},jj];
            end
        end
    end
    events{ii}(events_flags{ii},:) = [];
end

