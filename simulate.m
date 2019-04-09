

clear all;

niters = 10;

data = load_data;
load results_glme_fig3_nozscore.mat;

w = getEffects(results_VTURU, true);
perf_VTURU = gen(data, w, niters);

w = getEffects(results_V, true, [1 0 0]);
perf_V = gen(data, w, niters);

w = getEffects(results_VRU, true, [2 1 0]);
perf_VRU = gen(data, w, niters);

w = getEffects(results_VTU, true, [0 0 1]);
perf_VTU = gen(data, w, niters);

save simulate.mat
