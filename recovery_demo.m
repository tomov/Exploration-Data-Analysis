x = rand(100,1) - 0.5;
w = 3;
y = w * x;
c = binornd(1, normcdf(y));
tbl = table(c, x);
formula = 'c ~ -1 + x';

res = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal');

fprintf('simulated w = %.3f, fitted w = %.3f\n', w, fixedEffects(res));
