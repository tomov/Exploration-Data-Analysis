% toy example sanity checking univariate decoder
% => it works

clear;

n = 10; % # TRs
m = 5; % # regressors

X = rand(n,m) * 10;
b = rand(m,1) - 0.5;

y = X * b + rand(n,1);

bhat = glmfit(X,y, 'normal','constant', 'off')

reg = 3;

dec = (y - X(:,[1:reg-1 reg+1:end]) * bhat([1:reg-1 reg+1:end],:)) / bhat(reg);

[r,p] = corr(dec, X(:,reg));

r
p
