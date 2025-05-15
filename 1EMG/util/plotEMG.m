function plotEMG(traces, time, bound1, bound2)
num_traces = size(traces,1);

if (time == 0)
    hold on
    for i=1:num_traces
        plot(traces(i,:), '-k');
    end
    hold off;
else
    hold on
    for i=1:num_traces
        plot(traces(i,time-bound1:time+bound2), '-k');
    end
    hold off;
end
end