data = load_data;

% no standardize (no z-score)

tbl = data2table(data,0,1); % don't standardize, exclude timeouts
tbl = tbl(mod(tbl.run, 2) == 1,:); % odd runs only

formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')

save results_glme_fig3_odd_nozscore results_VTURU


% standardize (z-score)

tbl = data2table(data,1,1); % standardize, exclude timeouts
tbl = tbl(mod(tbl.run, 2) == 1,:); % odd runs only

formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')

save results_glme_fig3_odd results_VTURU
