
function [results_VTURU] =  getModeAll(data, subj)

    tbl = data2table(data,1);

    
    formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
    results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace');
    

end 

    
    
