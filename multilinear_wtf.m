data = load_data;
tbl = data2table(data, 0, 1);


[r, p] = corr(tbl.RU, tbl.TU)


[r, p] = corr(abs(tbl.RU), tbl.TU)

%{
p = binocdf(sum(tbl.C), length(tbl.C), 0.5);
p


r1 = [];
r2 = [];
for i = 1:length(data)
    r1 = [r1; data(i).r1];
    r2 = [r2; data(i).r2];
end
[h, p, ci, stats] = ttest2(r1, r2)


formula_orig = 'C ~ -1 + V + RU + VTU';
results_orig = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')


formula_intercept = 'C ~ 1 + V + RU + VTU';
results_intercept = fitglme(tbl,formula_intercept,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
c = compare(results_orig, results_intercept);
c


signRU = (tbl.RU >= 0) * 1 + (tbl.RU < 0) * (-1);
tbl = [tbl table(signRU)];
formula_sign = 'C ~ -1 + V + RU + VTU + signRU';
results_sign = fitglme(tbl,formula_sign,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
c = compare(results_orig, results_sign);
c


trial = [];
for i = 1:length(data)
    trial = [trial; data(i).trial(~data(i).timeout)];
end
tbl = [tbl table(trial)];
formula_trial = 'C ~ -1 + V + RU + VTU + trial';
results_trial = fitglme(tbl,formula_trial,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
c = compare(results_orig, results_trial);
c
%}


load multilinear_analysis_RU_glmbadre_ridge_CV.mat

plot([tbl.RU(1:1000), decRU(1:1000)]); legend({'RU', 'decRU'});

[w, names, stats] = fixedEffects(results_both{c});
DV_orig = tbl.RU * w(1) + tbl.VTU * w(2) + tbl.V * w(3);
[w, names, stats] = fixedEffects(results_both{c});
DV_both = tbl.RU * w(1) + tbl.VTU * w(2) + tbl.V * w(3) + tbl.decRU * w(4);

plot([DV_orig(1:1000), DV_both(1:1000)]); legend({'DV_{orig}', 'DV_{both}'});
