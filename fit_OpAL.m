function [results] = fit_OpAL(data, nstarts, outfile, hierarchical, fixedEffects)

if nargin < 4
    hierarchical = 0;
end
if nargin < 5
    fixedEffects = 0;
end

% create parameter structure using weakly informative priors

param(1).name = 'G_0R';
param(1).logpdf = @(x) sum(log(gampdf(x,1,5)));  % log density function for prior
param(1).lb = 0;    % lower bound
param(1).ub = 50;   % upper bound

param(2).name = 'N_0R';
param(2).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(2).lb = 0;  
param(2).ub = 50;
  
param(3).name = 'G_0S';
param(3).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(3).lb = 0;
param(3).ub = 50;

param(4).name = 'N_0S';
param(4).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(4).lb = 0;  
param(4).ub = 50;

param(5).name = 'V_0';
param(5).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(5).lb = 0;    
param(5).ub = 50; 

param(6).name = 'alpha';
param(6).logpdf = @(x) sum(log(betapdf(x,1.2,1.2)));
param(6).lb = 0;
param(6).ub = 1;

param(7).name = 'a';
param(7).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(7).lb = 0;
param(7).ub = 50;

param(8).name = 'b';
param(8).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(8).lb = 0;
param(8).ub = 50; 

for s = 1:length(data)
    data(s).N = sum(~data(s).timeout); % we ignore timeouts when computing loglik
end

if fixedEffects
    assert(~hierarchical);

    new_data.block = [];
    new_data.choice = [];
    new_data.reward = [];
    new_data.timeout = [];
    new_data.cond = [];
    new_data.N = 0;
    for s = 1:length(data)
        new_data.block = [new_data.block; data(s).block];
        new_data.reward = [new_data.reward; data(s).reward];
        new_data.choice = [new_data.choice; data(s).choice];
        new_data.timeout = [new_data.timeout; data(s).timeout];
        new_data.cond = [new_data.cond; data(s).cond];
        new_data.N = new_data.N + data(s).N;
    end
    data = new_data;
    save('shit.mat');
end


if hierarchical
    results = mfit_optimize_hierarchical(@loglik_OpAL, param, data, nstarts);
else
    results = mfit_optimize(@loglik_OpAL, param, data, nstarts);
end

save(outfile, 'results', 'data', 'param', 'nstarts', 'hierarchical', 'fixedEffects');
