function stacked = bin_and_stack(trains, duration)

nROIs = length(trains);
nUnits = [length(trains{1}), length(trains{2})];

binTrains = {};
for t=1:nROIs
    binTrains{t} = [];
    for i=1:nUnits(t)
        this_train = trains{t}{i};
        bt = train_to_binary(this_train,ceil(duration));
        binTrains{t} = [binTrains{t} ; bt]; %just concatenating 
    end
end

stacked = [binTrains{1} ; binTrains{2}]; %stacking!

disp('done converting to binary trains!!');

