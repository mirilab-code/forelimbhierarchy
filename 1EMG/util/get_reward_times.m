function rw = get_reward_times(digital, reward_ind)

ind1 = reward_ind(1);
ind2 = reward_ind(2);
ind3 = reward_ind(3);
ind4 = reward_ind(4);


d1 = diff(digital(ind1,:));
d2 = diff(digital(ind2,:));
d3 = diff(digital(ind3,:));
d4 = diff(digital(ind4,:));

rw1 = find(d1>0);
rw2 = find(d2>0);
rw3 = find(d3>0);
rw4 = find(d4>0);

rw = cell(4,1);
rw{ind1} = rw1;
rw{ind2} = rw2;
rw{ind3} = rw3;
rw{ind4} = rw4;

end