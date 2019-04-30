function [pred, mse] = multilinear_fit(X, y, Xtest, method, foldid, exclude)

    if exist('exclude', 'var')
        X = X(~exclude, :);
        y = y(~exclude);
        Xtest = Xtest(~exclude, :);
        foldid = foldid(~exclude);
    end

    % TODO all MSE's should be computed with CV (like fitnet)
    % TODO rm Xtest -- sometimes we predict based on X only ... (e.g. CV_CV)

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
            Lambda = 1e5;
            coef = ridge(y, X, Lambda, 0);
            pred = [ones(size(Xtest, 1), 1), Xtest] * coef; % include intercept term
            mse = immse(y, [ones(size(X, 1), 1), X] * coef);

        case 'ridge_CV'
            cv = cvpartition_from_folds(foldid);

            % find best lambda using CV
            Lambda = logspace(-10,20,30);
            m = [];
            for i = 1:length(Lambda)
                f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                m(i) = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
            end
            [~, idx] = min(m);
            fprintf('      min lambda(%d) = %f\n', idx, Lambda(idx));

            % compute predictions
            pred = ridgepred(X, y, Xtest, Lambda(idx));

            % use CV MSE...
            mse = m(idx); %immse(y, ridgepred(X, y, X, Lambda(idx))); <-- ...NOT the whole-data MSE

        case 'ridge_CV_CV'
            % actually use CV predictions & results
            %
            cv = cvpartition_from_folds(foldid);
            Lambda = logspace(-10,10,21);
            m = [];
            for i = 1:length(Lambda)
                f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                m(i) = crossval('mse', X, y, 'Predfun', f);
            end
            [~, idx] = min(m);
            %fprintf('                                                                  min lambda = %d (%e)\n', idx, Lambda(idx));

            f = @(XTRAIN,ytrain,XTEST,ytest) ridgepred(XTRAIN,ytrain,XTEST, Lambda(idx));

            % TODO FIXME this is wrong -- it expects you to compute e.g. a classification accuracy for each fold, not the actual predictions.. see docs
            % see elasticnet_CV_CV
            pred = crossval(f, X, y);
            pred = pred';
            pred = pred(:);

            mse = immse(y, pred);
            %mse = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
            %fprintf('                                                                  mse sanity: %e vs. %e\n', mse, immse(pred, y));


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

        case 'fitnet'

            % compute predictions
            pred = fitnet_pred(X, y, Xtest);

            % compute CV MSE
            cv = cvpartition_from_folds(foldid);
            f = @(XTRAIN,ytrain,XTEST) fitnet_pred(XTRAIN,ytrain,XTEST);
            mse = crossval('mse', X, y, 'Predfun', f);

        case 'elasticnet_CV'

            % hybrid ridge / lasso regression
            % lambda picked using CV
            cv = cvpartition_from_folds(foldid);
            [B, FitInfo] = lasso(X, y, 'Alpha', 0.5, 'CV', cv, 'NumLambda', 20);

            % no CV for prediciton
            idx = FitInfo.Index1SE;
            coef = B(:, idx);
            coef0 = FitInfo.Intercept(idx);

            pred = Xtest * coef + coef0;

            save wtf.mat
            mse = immse(y, X * coef + coef0);

        case 'elasticnet_CV_CV'

            % CV both to pick lambda, and to generate predictions later

            % pick lambda
            cv = cvpartition_from_folds(foldid);
            [B, FitInfo] = lasso(X, y, 'Alpha', 0.5, 'CV', cv, 'NumLambda', 20);
            Lambda = FitInfo.Lambda1SE;

            % predict
            kfold = length(unique(foldid));
            for k = 1:kfold
                which = training(cv, k);
                [B, FitInfo] = lasso(X(which,:), y(which,:), 'Alpha', 0.5, 'Lambda', Lambda);

                coef = B(:, 1);
                coef0 = FitInfo.Intercept(1);

                pred(~which,:) = X(~which,:) * coef + coef0;
            end

            mse = immse(pred, y);

        otherwise
            assert(false);
    end

end



function pred = ridgepred(X, y, Xtest, Lambda)
    coef = ridge(y, X, Lambda, 0);
    Xtest = [ones(size(Xtest, 1), 1), Xtest]; % include intercept term
    pred = Xtest * coef;
end

function pred = fitnet_pred(X, y, Xtest)
    X = X';
    y = y';
    Xtest = Xtest';
    net = fitnet(10);
    net = configure(net, X, y);
    net.trainParam.showWindow = false;
    net = train(net, X, y);

    pred = net(Xtest);
    pred = pred';
end

