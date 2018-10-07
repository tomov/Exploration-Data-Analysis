function [pred, mse] = multilinear_fit(X, y, Xtest, method, foldid)

    % Fit y = X * b using different methods.
    % Return mean-squared error (mse) and pred = Xtest * b.
    % foldid is only required for cross-validation
    %
    % USAGE:
    %   [pred, mse] = multilinear_fit(X, y, Xtest, method, foldid)
    %

    switch method
        case 'fitlm'
            mdl = fitlm(X, y, 'Intercept', true);
            pred = predict(mdl, Xtest);
            mse = mdl.MSE;

        case 'fitrlinear_ridge'
            mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge');
            mdl
            pred = predict(mdl, Xtest);
            mse = loss(mdl, X, y);

        case 'fitrlinear_ridge_CV'
            % find good lambda using CV
            cv = cvpartition_from_folds(foldid);
            Lambda = logspace(-10,10,21);
            cvmdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'CVPartition', cv, 'Lambda', Lambda);
            [l, idx] = min(kfoldLoss(cvmdl));
            %fprintf('      min lambda = %d\n', idx);

            mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', Lambda(idx));
            pred = predict(mdl, Xtest);
            mse = loss(mdl, X, y);

        case 'fitrlinear_lasso'
            mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'lasso');
            pred = predict(mdl, Xtest);
            mse = loss(mdl, X, y);

        case 'fitrlinear_lasso_CV'
            % find good lambda using CV
            cv = cvpartition_from_folds(foldid);
            Lambda = logspace(-10,10,21);
            cvmdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'lasso', 'CVPartition', cv, 'Lambda', Lambda);
            [l, idx] = min(kfoldLoss(cvmdl));
            %fprintf('      min lambda = %d\n', idx);

            mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'lasso', 'Lambda', Lambda(idx));
            pred = predict(mdl, Xtest);
            mse = loss(mdl, X, y);

        case 'ridge'
            Lambda = 0.1;
            coef = ridge(y, X, Lambda, 0);
            pred = [ones(size(Xtest, 1), 1), Xtest] * coef; % include intercept term
            mse = immse(y, [ones(size(X, 1), 1), X] * coef);

        case 'ridge_CV'
            cv = cvpartition_from_folds(foldid);
            Lambda = logspace(-10,10,21);
            m = [];
            for i = 1:length(Lambda)
                f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                m(i) = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
            end
            [~, idx] = min(m);
            %fprintf('      min lambda = %d\n', idx);

            pred = ridgepred(X, y, Xtest, Lambda(idx));
            mse = immse(y, ridgepred(X, y, X, Lambda(idx)));

        case 'lasso'
            [coef, FitInfo] = lasso(X, y, 'NumLambda', 1);
            coef0 = FitInfo.Intercept;
            pred = Xtest * coef + coef0;
            mse = immse(y, X * coef + coef0);

        case 'lasso_CV'
            cv = cvpartition_from_folds(foldid);
            [B, FitInfo] = lasso(X, y, 'CV', cv);
            idx = FitInfo.Index1SE;
            coef = B(:, idx);
            coef0 = FitInfo.Intercept(idx);
            pred = Xtest * coef + coef0;
            mse = immse(y, X * coef + coef0);

        otherwise
            assert(false);
    end

end



function pred = ridgepred(X, y, Xtest, Lambda)
    coef = ridge(y, X, Lambda, 0);
    Xtest = [ones(size(Xtest, 1), 1), Xtest]; % include intercept term
    pred = Xtest * coef;
end
