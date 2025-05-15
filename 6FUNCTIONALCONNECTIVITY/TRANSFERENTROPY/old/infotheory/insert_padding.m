function newtrain = insert_padding(train,times,pad)

padding = zeros(1,pad);
newtrain = [];
times = [1 times length(times)];
for i=2:length(times)
    if(i==2)
        a = train(times(i-1):times(i));
    else 
        a = train(times(i-1)+1:times(i));
    end

    newtrain = [newtrain a padding];
    
end

last = train( times(length(times)-1)+1:end);
newtrain = [newtrain last];

end