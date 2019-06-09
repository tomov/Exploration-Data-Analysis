function [results_V, results_VTU, results_VRU, results_VTURU ] = model_comparison_random(data)

    % same as model_comparison but with random starts like univariate_decoder_refactored
    % => also try a bunch b/c it fails, and take bestof
    % notice that here we DO INCLUDE bad_runs, and we do want that b/c the more behavior the better
    %
    
    tbl = data2table(data,0,1); % don't standardize, exclude timeouts

    best_of = 3; % get best model (BIC-wise) out of how many
    
    formula = 'C ~ -1 + V + (-1 + V|S)';
    successes = 0;
    for attempt = 1:100
        try
            res = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal', 'StartMethod', 'random')

            if ~exist('results_V', 'var') || results_V.LogLikelihood < res.LogLikelihood
                results_V = res;
            end
            successes = successes + 1;

            if successes == best_of
                break;
            end
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times...');



    formula = 'C ~ -1 + VTU + (-1 + VTU|S)';
    successes = 0;
    for attempt = 1:100
        try
            res = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal', 'StartMethod', 'random')

            if ~exist('results_VTU', 'var') || results_VTU.LogLikelihood < res.LogLikelihood
                results_VTU = res;
            end
            successes = successes + 1;

            if successes == best_of
                break;
            end
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times...');


    
    formula = 'C ~ -1 + V + RU + (-1 + V + RU|S)';
    successes = 0;
    for attempt = 1:100
        try
            res = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal', 'StartMethod', 'random')

            if ~exist('results_VRU', 'var') || results_VRU.LogLikelihood < res.LogLikelihood
                results_VRU = res;
            end
            successes = successes + 1;

            if successes == best_of
                break;
            end
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times...');
   



    formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
    successes = 0;
    for attempt = 1:100
        try
            res = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal', 'StartMethod', 'random')

            if ~exist('results_VTURU', 'var') || results_VTURU.LogLikelihood < res.LogLikelihood
                results_VTURU = res;
            end
            successes = successes + 1;

            if successes == best_of
                break;
            end
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times...');



    %save results_glme_fig3_nozscore_TrustRegion2D results_V results_VTU results_VRU results_VTURU
    save results_glme_fig3_nozscore_TrustRegion2D_random results_V results_VTU results_VRU results_VTURU
