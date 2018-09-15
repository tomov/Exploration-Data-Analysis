function [results] = fit_OpAL(data, nstarts, outfile, hierarchical, fixedEffects)

if nargin < 4
    hierarchical = 0;
end
if nargin < 5
    fixedEffects = 0;
end

% create parameter structure using weakly informative priors

param(1).name = 'alpha';
param(1).logpdf = @(x) sum(log(betapdf(x,1.2,1.2)));
param(1).lb = 0;
param(1).ub = 1;

param(2).name = 'a';
param(2).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(2).lb = 0;
param(2).ub = 50;

param(3).name = 'b';
param(3).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(3).lb = 0;
param(3).ub = 50; 

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
end


switch hierarchical
    case 2
        results = hfit_optimize(@loglik_OpAL, param, data, nstarts);
    case 1
        results = mfit_optimize_hierarchical(@loglik_OpAL, param, data, nstarts);
    otherwise
        results = mfit_optimize(@loglik_OpAL, param, data, nstarts);
end

save(outfile, 'results', 'data', 'param', 'nstarts', 'hierarchical', 'fixedEffects');
