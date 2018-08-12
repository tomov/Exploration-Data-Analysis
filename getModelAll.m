function [results_VTURU] =  getModeAll(data, subj)

    tbl = data2table(data,1);

    disp(subj);
    
    formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
    results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'PLIterations',200, 'CovariancePattern','diagonal');

end 

    
    
