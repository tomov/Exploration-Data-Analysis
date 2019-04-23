% Compare log likelihood computed manually vs. from GLM
% => sometimes they differ!!!!!!!!!! omg...
%

X = rand(10,1);
y = binornd(1, normcdf(X));
S = [1 1 1 1 1 2 2 2 2 2]';

T = table(y,X,S);

res = fitglme(T, 'y ~ X + (X | S)', 'Distribution', 'Binomial', 'Link', 'Probit', 'FitMethod', 'Laplace');

yhat = predict(res, T, 'Conditional', false); % <-- from fixed effects only
yhat2 = predict(res, T, 'Conditional', true); % <-- from fixed & random effects

lik = binopdf(y, 1, yhat);
lik2 = binopdf(y, 1, yhat2);

fprintf('Loglik from fixed effects: %.4f\n', sum(log(lik)));
fprintf('Loglik from fixed & random effects: %.4f\n', sum(log(lik2)));
fprintf('Loglik from model: %.4f\n', res.LogLikelihood);


