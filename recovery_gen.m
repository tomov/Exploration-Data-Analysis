% same as recovery.m but actually generative (i.e. don't take into account subject choices but generate our own)

clear all;

niters = 1000;

data = load_data;

formula = 'C ~ -1 + V + RU + VTU'; % single subject
w_orig = [];
w_rec = [];

for iter = 1:niters
    w = mvnrnd([0 0 0], 1 * eye(3));
    %w = exprnd([10 10 10]);
    disp(w);

    tbl = data2table_gen(data,0,1,w); % run generatively

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

save recovery_gen_mvnrnd.mat
