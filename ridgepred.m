% train & test ridge regression on single fold
%
function pred = ridgepred(X, y, Xtest, Lambda)
    coef = ridge(y, X, Lambda, 0);
    Xtest = [ones(size(Xtest, 1), 1), Xtest]; % include intercept term
    pred = Xtest * coef;
end

