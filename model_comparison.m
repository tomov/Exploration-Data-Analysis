function [results_V, results_VTU, results_VRU, results_VTURU ] = model_comparison(data)
    
    tbl = data2table(data,0,1); % don't standardize, exclude timeouts
    
    formula = 'C ~ -1 + V + (-1 + V|S)';
    results_V = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + VTU + (-1 + VTU|S)';
    results_VTU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + (-1 + V + RU|S)';
    results_VRU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
    results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'EBMethod', 'TrustRegion2D', 'CovariancePattern','diagonal')
    
    %save results_glme_fig3_nozscore_TrustRegion2D results_V results_VTU results_VRU results_VTURU
