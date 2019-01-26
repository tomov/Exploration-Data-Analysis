% test parameter recovery capabilities of model
% see https://psyarxiv.com/46mbn/

clear all;

niters = 1000;

data = load_data;

tbl = data2table(data,1,1); % standardize, exclude timeouts (behavior)
V = tbl.V;
RU = tbl.RU;
VTU = tbl.VTU;

formula = 'C ~ -1 + V + RU + VTU'; % single subject
w_orig = [];
w_rec = [];

for iter = 1:niters
    w = mvnrnd([0 0 0], 10 * eye(3));
    %w = exprnd([10 10 10]);
    disp(w);

    DV = w(1) * V + w(2) * RU + w(3) * VTU;
    pred = normcdf(DV); % manual prediction

    C = binornd(1, pred);
    tbl.C = C;

    try
        results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal');

        w_orig = [w_orig; w];
        [w, names] = fixedEffects(results_VTURU);
        w_rec = [w_rec; w(3) w(1) w(2)];
    catch e
        disp('got an error while fitting...');
        disp(e);
        % TODO might introduce correlations between parameters
    end
end

save recovery_mvnrnd.mat
