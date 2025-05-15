function si = spike_info(isi,pd)
%     minnotzero = min(isi(isi > 0));
    si = [];
    for i=1:length(isi)
        this_isi = isi(i);
%         if(this_isi == 0)
%             this_isi = minnotzero;
%         end
        this_info = pdf(pd,this_isi);

        si = [si this_info];
    end
end