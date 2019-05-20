function [pred, mse, mses] = multilinear_fit(X, y, Xtest, method, foldid, exclude, Lambda)

    if exist('exclude', 'var') && ~isempty(exclude)
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

    % remap foldid's from 1..k
    %
    fid = unique(foldid);
    for k = 1:length(fid)
        foldid(foldid == fid(k)) = k;
    end

    mses = [];
    pred = [];
    mse = [];

    switch method
        case 'fitlm'
            mdl = fitlm(X, y, 'Intercept', true);
            pred = predict(mdl, Xtest);
            mse = mdl.MSE;

        case 'fitrlinear_ridge' % better than ridge for high-dimensional data
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

        case 'fitrlinear_CV_1' % like ridge_CV_1

            % compute MSE for different lambdas

            cv = cvpartition_from_folds(foldid);

            if ~exist('Lambda', 'var')
                Lambda = logspace(-10,20,30);
            end

            cvmdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'CVPartition', cv, 'Lambda', Lambda);

            mses = kfoldLoss(cvmdl);


        case 'fitrlinear_CV_2' % like ridge_CV_2

            % use given lambda for prediction
            % TODO redundant to re-fit -- we already fitted this model in fitrlinear_CV_1; just reuse it

            cv = cvpartition_from_folds(foldid);

            cvmdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'CVPartition', cv, 'Lambda', Lambda);

            pred = kfoldPredict(cvmdl);
            mse = kfoldLoss(cvmdl);

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

            % use CV to pick lambda & then generate predictions

            % pick lambda
            cv = cvpartition_from_folds(foldid);
            Lambda = logspace(-10,20,30);
            m = [];
            for i = 1:length(Lambda)
                f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                m(i) = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
            end
            [~, idx] = min(m);
            fprintf('                                                                  min lambda(%d) = %f\n', idx, Lambda(idx));

            % predict
            kfold = length(unique(foldid));
            for k = 1:kfold
                which = training(cv, k);
                pred(~which,:) = ridgepred(X(which,:), y(which,:), X(~which,:), Lambda(idx));
            end
            mse = immse(pred, y);

            fprintf('                                                                  mse sanity: %f vs. %f\n', mse, m(idx));
            if abs(mse - m(idx)) > 1e-10
                save fuck.mat
                assert(abs(mse - m(idx)) < 1e-10);
            end

        case 'ridge_CV_1'

            % first half of ridge_CV_CV -- try different lambdas

            cv = cvpartition_from_folds(foldid);

            if ~exist('Lambda', 'var')
                Lambda = logspace(-10,20,30);
            end

            for i = 1:length(Lambda)
                f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                mses(i) = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
            end

        case 'ridge_CV_2'

            % second half of ridge_CV_CV -- use given lambda to predict

            cv = cvpartition_from_folds(foldid);

            kfold = length(unique(foldid));
            for k = 1:kfold
                which = training(cv, k);
                pred(~which,:) = ridgepred(X(which,:), y(which,:), X(~which,:), Lambda);
            end
            mse = immse(pred, y);


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

            % TODO it's broken for subject 3, FitInfo has NaNs, so can't select lambda .... w t f

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

            % TODO it's broken for subject 3, FitInfo has NaNs, so can't select lambda .... w t f

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

