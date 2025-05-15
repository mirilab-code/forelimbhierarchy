function mX = get_trials_and_avg(M,times,window)

n = size(M,1);
X = zeros(n,length(window),length(times));
exc = [];
for i=1:length(times)
    t = times(i);
    try
        chunk = M(:,t+window);
        X(:,:,i) = chunk;
    catch
        exc = [exc i];
    end
end

X(:,:,exc) = [];

mX = mean(X,3);

end