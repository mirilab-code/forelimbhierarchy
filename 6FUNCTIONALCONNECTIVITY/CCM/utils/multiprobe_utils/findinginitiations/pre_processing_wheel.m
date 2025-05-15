function [climbing_struct,filtered_wheel_signal,isclimbing]=...
    pre_processing_wheel(raw_wheel_signal,immobility_min_duration,immobility_min_duration_before_climbing,min_climbing_angle)

sr = 1000;      % at least for all the mouse palooza stuff
% sr: sampling rate in Hz
% immobility_min_duration %minimal duration of "immobility" periods (sec)
% immobility_min_duration_before_climbing %minimal duration of the immobility period preceding "climbing" bouts (sec)
% a "mobility" fragment which is not passing this criterion will not be considered as "climbing"
% min_climbing_angle %minimum angular distance travelled during a climbing bout (°)

%% filter signal
pos_pics=find(diff(raw_wheel_signal)>3);
fragment_index(1,1)=1;
fragment_index(1,2)=pos_pics(1);
for i_fragment=1:numel(pos_pics)-1
    fragment_index(i_fragment+1,1)=pos_pics(i_fragment)+1;
    fragment_index(i_fragment+1,2)=pos_pics(i_fragment+1);
end
fragment_index(size(fragment_index,1)+1,:)=[pos_pics(end)+1,numel(raw_wheel_signal)]; %fragment déterminé par saut de 0 à 3.31V
for i_fragment=1:size(fragment_index,1)
    filtered_wheel_signal{i_fragment,1}=  smoothdata(raw_wheel_signal(fragment_index(i_fragment,1):fragment_index(i_fragment,2)),'sgolay',1000,'Degree',1);
end
filtered_wheel_signal=cell2mat(filtered_wheel_signal');

%% Detection immobility periods, duration of immobility periods, mobility periods & travelled angle

climbed_angle=generate_climbed_angle(filtered_wheel_signal,sr,max(filtered_wheel_signal));
angular_derivative=climbed_angle(:,2);

size_window=round(.1*sr); %ms

for_movement_detection=movsum(angular_derivative,[0 size_window],"Endpoints","discard");

% figure
% yyaxis left
% plot(raw_wheel_signal)
% yyaxis right
% plot(abs(for_movement_detection))

immobile=abs(for_movement_detection)<.5; %°
non_immobile=logical(1-immobile);

% figure
% yyaxis left
% plot((1:numel(immobile))/sr,filtered_wheel_signal(size_window+1:numel(angular_derivative)))
% yyaxis right
% plot((1:numel(immobile))/sr,immobile)

indices_periodes_immobilite=find(immobile)';
indices_periodes_mobilite=find(non_immobile)';
indices_periodes_immobilite2=mat2cell(indices_periodes_immobilite,1,diff([0,find(diff(indices_periodes_immobilite)~=1),length(indices_periodes_immobilite)]));
duree_periodes_immobilite=cellfun(@numel,indices_periodes_immobilite2)./sr;
% supprimer périodes d'immobilité de moins de immobility_min_duration s et les considérer comme périodes mobiles:
indices_a_exclure_dimmobilite=cell2mat(indices_periodes_immobilite2(duree_periodes_immobilite<immobility_min_duration));
indices_periodes_mobilite=sort([indices_periodes_mobilite,indices_a_exclure_dimmobilite]);
indices_periodes_mobilite2=mat2cell(indices_periodes_mobilite,1,diff([0,find(diff(indices_periodes_mobilite)~=1),length(indices_periodes_mobilite)]));
indices_periodes_immobilite=setdiff(indices_periodes_immobilite,indices_a_exclure_dimmobilite);
indices_periodes_immobilite2=mat2cell(indices_periodes_immobilite,1,diff([0,find(diff(indices_periodes_immobilite)~=1),length(indices_periodes_immobilite)]));
% duree_periodes_immobilite=cellfun(@numel,indices_periodes_immobilite2)./sr;

angle_travelled=cellfun(@(x) sum(angular_derivative(x(1)+size_window+1:x(end)+size_window)),indices_periodes_mobilite2);
%supprimer périodes de mobilité dont distance cumulée est inférieure à min_climbing_angle° mais ne pas les considérer comme périodes immobiles:
indices_periodes_mobilite2(angle_travelled<min_climbing_angle)=[];
angle_travelled=cellfun(@(x) sum(angular_derivative(x(1)+size_window+1:x(end)+size_window)),indices_periodes_mobilite2);

indices_periodes_mobilite2=cellfun(@(x) x+size_window,indices_periodes_mobilite2,'UniformOutput',false);
indices_periodes_immobilite2=cellfun(@(x) x+size_window,indices_periodes_immobilite2,'UniformOutput',false);
indices_periodes_immobilite3=horzcat(indices_periodes_immobilite2{:});

% disregard mobility periods which are not preceded by immobility lasting at least immobility_min_duration_before_climbing s:
if immobility_min_duration_before_climbing~=0
    to_be_deleted=[];
    for i_period=1:numel(indices_periodes_mobilite2)
        if indices_periodes_mobilite2{i_period}(1)-round(immobility_min_duration_before_climbing*sr)>=1
            frames_during_which_should_be_immobile=...
                indices_periodes_mobilite2{i_period}(1)-round(immobility_min_duration_before_climbing*sr):indices_periodes_mobilite2{i_period}(1)-round(.01*sr);
        else
            frames_during_which_should_be_immobile=...
               1:indices_periodes_mobilite2{i_period}(1)-round(.01*sr);
        end
        if ~all(ismember(frames_during_which_should_be_immobile,indices_periodes_immobilite3))
            to_be_deleted=[to_be_deleted,i_period];
        end
    end
end
climbing=indices_periodes_mobilite2';
climbing(to_be_deleted)=[];
angle_travelled(to_be_deleted)=[];



climbing_onsets = cellfun(@(x) x(1), climbing);
climbing_offsets = cellfun(@(x) x(end), climbing);
climbing_durations = cellfun(@(x) length(x), climbing);

climbing_struct = struct('epoch_times',climbing);
onsets = num2cell(climbing_onsets);
[climbing_struct(:).onset] = deal(onsets{:});
offsets = num2cell(climbing_offsets);
[climbing_struct(:).offset] = deal(offsets{:});
durs = num2cell(climbing_durations);
[climbing_struct(:).duration] = deal(durs{:});
at = num2cell(angle_travelled);
[climbing_struct(:).angle_traveled] = deal(at{:});


%% plots

climbing_for_plot=false(1,numel(filtered_wheel_signal));
for i_period=1:numel(indices_periodes_mobilite2)
    climbing_for_plot(indices_periodes_mobilite2{i_period})=true;
end

isclimbing = climbing_for_plot;

x = 1:length(filtered_wheel_signal);
figure
hold on;
plot(x,filtered_wheel_signal)
plot(x(climbing_for_plot),filtered_wheel_signal(climbing_for_plot),'.r')
hold off

% for i_period=1:numel(indices_periodes_mobilite2)
%     angular_derivative_climbing(:,i_period)=...
%         angular_derivative(indices_periodes_mobilite2{i_period}(1)-round(1*sr):indices_periodes_mobilite2{i_period}(1)+round(1*sr));
% end
% 
% fig_angularderiv = figure;
% plot(median(angular_derivative_climbing,2))
% title('median angular derivative')

end