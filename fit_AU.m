function [results] = fit_AU(data, nstarts, outfile, hierarchical, fixedEffects)

if nargin < 4
    hierarchical = 0;
end
if nargin < 5
    fixedEffects = 0;
end

% create hyperparameter structure

hyparam(1).name = 'G_0R hyperparameters';
hyparam(1).lb = [0 0];
hyparam(1).ub = [10 10];
hyparam(1).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(1).rnd = @() [rand * 10 rand * 10];

hyparam(2).name = 'N_0R hyperparameters';
hyparam(2).lb = [0 0];
hyparam(2).ub = [10 10];
hyparam(2).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(2).rnd = @() [rand * 10 rand * 10];

hyparam(3).name = 'G_0S hyperparameters';
hyparam(3).lb = [0 0];
hyparam(3).ub = [10 10];
hyparam(3).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(3).rnd = @() [rand * 10 rand * 10];

hyparam(4).name = 'N_0S hyperparameters';
hyparam(4).lb = [0 0];
hyparam(4).ub = [10 10];
hyparam(4).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(4).rnd = @() [rand * 10 rand * 10];

hyparam(5).name = 'alpha hyperparameters';
hyparam(5).lb = [1 1];
hyparam(5).ub = [10 10];
hyparam(5).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(5).rnd = @() [rand * 10 rand * 10];

hyparam(6).name = 'beta hyperparameters';
hyparam(6).lb = [1 1];
hyparam(6).ub = [10 10];
hyparam(6).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(6).rnd = @() [rand * 10 rand * 10];

hyparam(7).name = 'a hyperparameters';
hyparam(7).lb = [0 0];
hyparam(7).ub = [10 10];
hyparam(7).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(7).rnd = @() [rand * 10 rand * 10];

hyparam(8).name = 'b hyperparameters';
hyparam(8).lb = [0 0];
hyparam(8).ub = [10 10];
hyparam(8).logpdf = @(x) log(unifpdf(x(1), 0, 10)) + log(unifpdf(x(2), 0, 10)); 
hyparam(8).rnd = @() [rand * 10 rand * 10];

% create parameter structure using weakly informative priors

param(1).name = 'G_0R';
param(1).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(1).hrnd = @(h) gamrnd(h(1),h(2));
param(1).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(1).lb = 0;
param(1).ub = 50;

param(2).name = 'N_0R';
param(2).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(2).hrnd = @(h) gamrnd(h(1),h(2));
param(2).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(2).lb = 0;  
param(2).ub = 50;
  
param(3).name = 'G_0S';
param(3).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(3).hrnd = @(h) gamrnd(h(1),h(2));
param(3).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(3).lb = 0;
param(3).ub = 50;

param(4).name = 'N_0S';
param(4).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(4).hrnd = @(h) gamrnd(h(1),h(2));
param(4).logpdf = @(x) sum(log(gampdf(x,1,5))); 
param(4).lb = 0;  
param(4).ub = 50;

param(5).name = 'alpha';
param(5).hlogpdf = @(x,h) sum(log(betapdf(x,h(1),h(2))));
param(5).hrnd = @(h) betarnd(h(1),h(2));
param(5).logpdf = @(x) sum(log(betapdf(x,1.2,1.2)));
param(5).lb = 0;
param(5).ub = 1;

param(6).name = 'beta';
param(6).hlogpdf = @(x,h) sum(log(betapdf(x,h(1),h(2))));
param(6).hrnd = @(h) betarnd(h(1),h(2));
param(6).logpdf = @(x) sum(log(betapdf(x,1.2,1.2)));
param(6).lb = 0;
param(6).ub = 1;

param(7).name = 'a';
param(7).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(7).hrnd = @(h) gamrnd(h(1),h(2));
param(7).logpdf = @(x) sum(log(gampdf(x,1,5)));
param(7).lb = 0;
param(7).ub = 50;

param(8).name = 'b';
param(8).hlogpdf = @(x,h) sum(log(gampdf(x,h(1),h(2))));
param(8).hrnd = @(h) gamrnd(h(1),h(2));
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

switch hierarchical
    case 2
        results = hfit_optimize(@loglik_AU, hyparam, param, data, nstarts);
    case 1
        results = mfit_optimize_hierarchical(@loglik_AU, param, data, nstarts);
    otherwise
        results = mfit_optimize(@loglik_AU, param, data, nstarts);
end

save(outfile, 'results', 'data', 'param', 'nstarts', 'hierarchical', 'fixedEffects');
