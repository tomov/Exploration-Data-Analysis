function perf = getPerfs(data) % get subject P(better option)

    perf = []; % P(better option)
    for s = 1:length(data)
        which = ~data(s).timeout;
        better = (data(s).mu2(which) > data(s).mu1(which)) + 1; 
        C = double(data(s).choice(which) == better); % human choices
        %perf = [perf; mean(data(s).reward)];
        perf = [perf; mean(C)];
    end
end
