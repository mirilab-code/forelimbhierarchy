function tsls = time_since_last_spike(isi)
    n = length(isi);
    tsls = [];
    for i=1:n
        this_isi = isi(i);
        chunk = 0:this_isi;
        tsls = [tsls chunk];
    end
end