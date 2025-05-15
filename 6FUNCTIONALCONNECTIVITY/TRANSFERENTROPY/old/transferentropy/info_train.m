function it = info_train(isi,info,pd)
    isi = reshape(isi,[],1);
    info = reshape(info,[],1);
    % an ISI longer than t is surprisingly long, an isi shorter is surprisingly short
    medianISI = find(info==min(info),1);
    medianInfo = info(medianISI);
    long = info(medianISI+1:end);

    nISI = length(isi);
    it = [];
    for i=1:nISI
        this_isi = floor(isi(i));
        if(this_isi == 0)
            this_isi = 1;
        end
        if(this_isi < medianISI)
            y = -log2(pdf(pd,this_isi));
%             disp(['short ' num2str(this_isi)]);
            chunk = zeros(this_isi,1)+medianInfo;
            chunk(end) = y;
%             disp([length(chunk) this_isi]);
        end
        if(this_isi >= medianISI)
%             disp(['long ' num2str(this_isi)]);
            chunk1 = zeros(medianISI,1)+medianInfo;
            chunk2 = long(1:this_isi-medianISI);
            chunk = [chunk1; chunk2];
%             disp([length(chunk) this_isi]);

        end
        
        it = [it; chunk];
    end

end