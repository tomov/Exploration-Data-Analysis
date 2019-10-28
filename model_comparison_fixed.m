function [results_V, results_VTU, results_VRU, results_VTURU, bics] = model_comparison_fixed(data)

    % same as model_comparison but fixed effects for each subject
    
    tbl_all = data2table(data,0,1); % don't standardize, exclude timeouts

    bics = [];

    for s = 1:length(data)

        bic = [];
        tbl = tbl_all(tbl_all.S == s, :);
    
        formula = 'C ~ -1 + V ';
        results_V{s} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal');
        bic = [bic results_V{s}.ModelCriterion.BIC];
        
        formula = 'C ~ -1 + VTU ';
        results_VTU{s} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal');
        bic = [bic results_VTU{s}.ModelCriterion.BIC];
        
        formula = 'C ~ -1 + V + RU ';
        results_VRU{s} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal');
        bic = [bic results_VRU{s}.ModelCriterion.BIC];
        
        formula = 'C ~ -1 + V + RU + VTU ';
        results_VTURU{s} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'EBMethod', 'TrustRegion2D', 'CovariancePattern','diagonal');
        bic = [bic results_VTURU{s}.ModelCriterion.BIC];
      
        %{
        formula = 'C ~ -1 + V + RU + VTU2 + (-1 + V + RU + VTU2|S)';
        results_VTU2RU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'EBMethod', 'TrustRegion2D', 'CovariancePattern','diagonal')

        formula = 'C ~ -1 + VTU2 + RUTU2 + (-1 + VTU2 + RUTU2|S)';
        results_VTU2RUTU2 = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'EBMethod', 'TrustRegion2D', 'CovariancePattern','diagonal')

        formula = 'C ~ -1 + VTU + RUTU + (-1 + VTU + RUTU|S)';
        results_VTURUTU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'EBMethod', 'TrustRegion2D', 'CovariancePattern','diagonal')

        %}

        bics = [bics; bic];
    end

    lme = -0.5 * bics;
    [alpha, exp_r, xp, pxp, bor] = bms(lme);

    pxp

    %save results_glme_fig3_nozscore_TrustRegion2D results_V results_VTU results_VRU results_VTURU
    save results_glme_fig3_nozscore_TrustRegion2D_fixed
