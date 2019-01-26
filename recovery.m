% test parameter recovery capabilities of model
% see https://psyarxiv.com/46mbn/

clear all;

niters = 100;

data = load_data;

tbl = data2table(data,1,1); % standardize, exclude timeouts (behavior)
V = tbl.V;
RU = tbl.RU;
VTU = tbl.VTU;

formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
w_orig = [];
w_rec = [];

for iter = 1:100
    w = mvnrnd([0 0 0], 10 * eye(3));

    DV = w(1) * V + w(2) * RU + w(3) * VTU;
    pred = normcdf(DV); % manual prediction

    C = pred > 0.5;
    tbl.C = C;

    w_orig = [w_orig; w];

    results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal');
    [w, names] = fixedEffects(results_VTURU);
    w_rec = [w_rec; w(3) w(1) w(2)];
end

save recovery.mat
