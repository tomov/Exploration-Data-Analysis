function simulate_more(nws, niters)


%nws = 16;
%niters = 10;

filename = sprintf('simulate_more_nws=%d_niters=%d.mat', nws, niters);
filename

data = load_data;
load results_glme_fig3_nozscore.mat;

w = getEffects(results_VTURU, true);

w1s = linspace(0, 1, nws);
w2s = linspace(0, 1, nws);
w3s = linspace(0, 1, nws);

m = nan(nws, nws, nws);
se = nan(nws, nws, nws);

for i = 1:nws
    i
    for j = 1:nws
        j
        for k = 1:nws
            w = [w1s(i) w2s(j) w3s(k)];
            perf = gen(data, w, niters);

            m(i,j,k) = mean(perf);
            se(i,j,k) = std(perf) / sqrt(length(perf));
        end
    end
end

save(filename, '-v7.3');
