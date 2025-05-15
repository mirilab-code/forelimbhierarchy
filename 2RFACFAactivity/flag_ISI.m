function bad_units =  flag_ISI(N)

thresh = 1;            % how small an ISI has to be to be considered an ISI violation
tolerance = 0.10;      % the percentage of ISI violations needed to flag the unit

trains = {N.train};
isi = cellfun(@diff, trains, 'UniformOutput', false);
flag = cellfun(@(trn) sum(trn<thresh)/length(trn),isi);

bad_units = find(flag>tolerance);











end