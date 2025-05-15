function spk_trn = create_fake_neuron(a,b,train_length)
    p = makedist('Gamma',a,b);
    
    % need to keep adding on isi's until the total time is <= duration
    spk_trn = [];
    while(sum(spk_trn) < train_length)
        this_isi = random(p,1);
        spk_trn = [spk_trn this_isi];
    end
end