function [results] = fit_OpAL(data, nstarts, outfile, hierarchical)

if nargin < 4
    hierarchical = 0;
end

% create parameter structure using weakly informative priors

param(1).name = 'G_0';
param(1).logpdf = @(x) sum(log(gampdf(x,1,5)));  % log density function for prior
param(1).lb = 0;    % lower bound
param(1).ub = 50;   % upper bound

param(2).name = 'N_0';
param(2).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(2).lb = 0;  
param(2).ub = 50;

param(3).name = 'V_0';
param(3).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(3).lb = 0;    
param(3).ub = 50; 

param(4).name = 'alpha';
param(4).logpdf = @(x) sum(log(betapdf(x,1.2,1.2)));
param(4).lb = 0;
param(4).ub = 1;

param(5).name = 'a';
param(5).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(5).lb = 0;
param(5).ub = 50;

param(6).name = 'b';
param(6).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(6).lb = 0;
param(6).ub = 50; 

for s = 1:length(data)
    data(s).N = sum(~data(s).timeout); % we ignore timeouts when computing loglik
end

if hierarchical
    results = mfit_optimize_hierarchical(@loglik_OpAL, param, data, nstarts);
else
    results = mfit_optimize(@loglik_OpAL, param, data, nstarts);
end

save(outfile, 'results', 'data', 'param', 'nstarts');
