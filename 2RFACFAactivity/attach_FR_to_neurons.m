function Nnew = attach_FR_to_neurons(N,FR)

Nnew = N;

units = [N.unit];
for i=1:length(units)
    fr = FR(i,:);
    if(sum(fr) ~= 0)
        Nnew(i).ifr = fr;
    end
end


end