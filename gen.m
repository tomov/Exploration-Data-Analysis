function perf = gen(data, w, niters)

    % run model generatively for niters, return perf

    perf = [];

    for iter = 1:niters
        tbl = data2table_gen(data, 0, 1, w);

        for s = 1:length(data)
            which = tbl.S == s;
            better = tbl.mu1(which) > tbl.mu2(which);
            C = tbl.C(which) == better;
            perf = [perf; mean(C)];
        end
    end
end
