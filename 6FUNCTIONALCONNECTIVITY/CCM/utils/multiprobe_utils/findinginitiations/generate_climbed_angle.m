function climbed_angle=generate_climbed_angle(wheel_trace_fragment,sr,max_voltage)

climbed_angle=...
    [repmat(1/sr,numel(wheel_trace_fragment)-1,1),...
    correction_wheel_angle(diff(   wheel_trace_fragment   ),max_voltage)*360/max_voltage];

end

function corrected_delta_voltage=correction_wheel_angle(delta_voltage,max_encoder_voltage)
% déduire distance angulaire parcourue, en considérant shifts voltage
% survenant lorsque roue passe point 3.6V=>0V.
% voltage diminue (jusqu'à 0) lorsque rouge avance (augmente lorsque reculte)
%             max_encoder_voltage couvre les 360° de rotation de la roue

corrected_delta_voltage=nan(numel(delta_voltage),1);
for i=1:numel(delta_voltage)
    if abs(delta_voltage(i))>=max_encoder_voltage-.6
        corrected_delta_voltage(i,1)=-sign(delta_voltage(i))*(abs(delta_voltage(i))-max_encoder_voltage);
    else
        corrected_delta_voltage(i,1)=-delta_voltage(i);
    end
end

end