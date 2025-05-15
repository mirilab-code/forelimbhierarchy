function [flagged_units,refract_pct] = checkViolations(events,acg_bin,max_time,probe)
%checkViolations returns flagged units with inter-spike interval
%violations. Violators are determined using autocorrelation (acf function).
%
% INPUTS:
%   events: events matrix from preprocessing
%   acg_bin: bin size in acf function (0.0001 s)
%   max_tme: max time in acf fucntion (0.1 s)
%   probe: probe number (only used for plotting)

flagged_units = [];
refract_pct = [];
train = events_to_train(events);

%Iterate through each neuron spike train
for ii = 1:length(train)
    if ~isempty(train{ii})
        disp(['Unit ' num2str(ii) ' of ' num2str(length(train))])
        %calculate the autocorrelation
        [ac,xbin] = acf(train{ii}/1000,acg_bin,max_time);
        ac(ceil(length(xbin)/2)) = NaN;
        
        %Define refractory period
        zero_index = ((length(xbin)+1)/2)+1;
        start_refract = zero_index + 3;
        end_refract = zero_index + 10;
        
        %Define flanking period
        start_flank = zero_index + 100;
        end_flank = zero_index + 500;
        
        %Sum values in these periods
        refract = sum(ac(start_refract:end_refract));
        flank = sum(ac(start_flank:end_flank));
        
        %Calculate refractory violation percentage as ratio of refractory
        %spikes to flanking spikes
        refract_per_sec = refract/(0.0001*(end_refract-start_refract+1));
        flank_per_sec = flank/(0.0001*(end_flank-start_flank+1));
        violation_pct = refract_per_sec/flank_per_sec;
        refract_pct = [refract_pct,violation_pct];
        
        % Check to see if violation percentage is above our threshold
        if violation_pct > .18
            flagged_units = [flagged_units,ii];
            
            % plotting ISI violator's acg
%             f = figure;
%             bar(xbin,ac,'r');
%             str = sprintf('Autocorrelogram: Probe %d: Unit %d',probe,ii);
%             title(str);
%             movegui(f,'south');
%             str = sprintf('Refract. counts/s: %5.1f\nEdges counts/s: %5.1f\nPercentage: %5.3f',...
%                 refract_per_sec,flank_per_sec,violation_pct);
%             annotation('textbox',[.15,.75,.35,.15],'String',str);
        else
            % plotting all other unit's acg
%             figure;
%             bar(xbin,ac,'b');
%             str = sprintf('Autocorrelogram: Probe %d: Unit %d',probe,ii);
%             title(str);
%             str = sprintf('Refract. counts/s: %5.1f\nEdges counts/s: %5.1f\nPercentage: %5.3f',...
%                 refract_per_sec,flank_per_sec,violation_pct);
%             annotation('textbox',[.15,.75,.35,.15],'String',str);
        end       
    end
end


