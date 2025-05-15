function [CR,RC] = separate(M,units)

RC = M(units(1)+1:end,1:units(1));
CR = M(1:units(1),units(1)+1:end);

end