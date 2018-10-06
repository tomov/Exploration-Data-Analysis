% multilinear regression analysis for RU for Badre 2012 RLPFC ROI
% try to decode |RU| from multivariate ROI activity and see if it predicts
% choices better than RU from model
%
% TODO dedupe with badre_2012_activations_analysis.m

function badre_2012_multilinear_analysis(method)

    printcode;

    EXPT = exploration_expt();

    data = load_data;

    formula_both = 'C ~ -1 + V + RU + VTU + decRU';
    formula_RU = 'C ~ -1 + V + RU + VTU';
    formula_decRU = 'C ~ -1 + V + decRU + VTU';

    filename = ['badre_2012_multilinear_analysis_', method, '.mat'];
    disp(filename);

    % clusters = masks from paper
    masks = badre_2012_create_masks(false);
    masks = masks(1); % TODO all masks

    % extract trial_onset (raw, unsmoothed) betas
    roi = extract_roi_betas(masks, 'trial_onset');
    save(filename, '-v7.3');

    load(filename, 'roi'); 

    [~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

    % clean up betas
    %
    for c = 1:length(roi)
        for s = 1:length(data)
            B = roi(c).subj(s).betas;
            runs = find(goodRuns{s});
            data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials
            which_nan = any(isnan(B(~data(s).exclude, :)), 1); % exclude nan voxels (ignoring bad runs and timeouts; we exclude those in the GLMs)
            B(:, which_nan) = [];
            data(s).betas{c} = B;
        end
    end

    % extract regressors
    %
    for s = 1:length(data)
        which_all = logical(ones(length(data(s).run), 1));
        [~, absRU] =  get_latents(data, s, which_all, 'abs');
        [~, RU] = get_latents(data, s, which_all, 'left');
        data(s).absRU = absRU;
        data(s).RU = RU; % for sign-correction
    end

    save(filename, '-v7.3');

    load(filename); 

    for c = 1:numel(masks)
        mask = masks{c};
        [~, masknames{c}, ~] = fileparts(mask);
        disp(mask);

        decRU = [];
        exclude = [];
        mse = [];
        for s = 1:length(data)
            exclude = [exclude; data(s).exclude];
            X = data(s).betas{c};
            y = data(s).absRU;

            % remove bad data points
            X = X(~data(s).exclude, :);
            y = y(~data(s).exclude);

            switch method
                case 'fitlm'
                    mdl = fitlm(X, y, 'Intercept', true);
                    pred = predict(mdl, data(s).betas{c}); % predict using full data set; we ignore bad trials later
                    mse(s) = mdl.MSE;

                case 'fitrlinear_ridge'
                    mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge');
                    mdl
                    pred = predict(mdl, data(s).betas{c}); % predict using full data set; we ignore bad trials later
                    mse(s) = loss(mdl, X, y);

                case 'fitrlinear_ridge_CV'
                    % find good lambda using CV
                    cv = cvpartition_from_folds(data(s).run(~data(s).exclude)); % one run per fold
                    Lambda = logspace(-5,0,100);
                    cvmdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'CVPartition', cv, 'Lambda', Lambda);
                    [l, idx] = min(kfoldLoss(cvmdl));
                    fprintf('      min lambda for %d = %f\n', s, idx);

                    mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', Lambda(idx));
                    pred = predict(mdl, data(s).betas{c}); % predict using full data set; we ignore bad trials later
                    mse(s) = loss(mdl, X, y);

                case 'fitrlinear_lasso'
                    mdl = fitrlinear(X, y, 'ObservationsIn', 'rows', 'Learner', 'leastsquares', 'Regularization', 'lasso');
                    pred = predict(mdl, data(s).betas{c}); % predict using full data set; we ignore bad trials later
                    mse(s) = loss(mdl, X, y);

                case 'ridge'
                    Lambda = 0.1;
                    coef = ridge(y, X, Lambda, 0);
                    pred = [ones(size(data(s).betas{c}, 1), 1), data(s).betas{c}] * coef; % include intercept term
                    mse(s) = immse(y, pred(~data(s).exclude));

                case 'ridge_CV'
                    cv = cvpartition_from_folds(data(s).run(~data(s).exclude)); % one run per fold
                    Lambda = logspace(-5,0,100);
                    m = [];
                    for i = 1:length(Lambda)
                        f = @(XTRAIN,ytrain,XTEST) ridgepred(XTRAIN,ytrain,XTEST, Lambda(i));
                        m(i) = crossval('mse', X, y, 'Predfun', f, 'partition', cv);
                    end
                    [~, idx] = min(m);
                    fprintf('      min lambda for %d = %f\n', s, idx);

                    pred = ridgepred(X, y, data(s).betas{c}, Lambda(idx));
                    mse(s) = immse(y, pred(~data(s).exclude));

                case 'lasso'
                    [coef, FitInfo] = lasso(X, y, 'NumLambda', 1);
                    coef0 = FitInfo.Intercept;
                    pred = data(s).betas{c} * coef + coef0;
                    mse(s) = immse(y, pred(~data(s).exclude));

                case 'lasso_CV'
                    cv = cvpartition_from_folds(data(s).run(~data(s).exclude)); % one run per fold
                    [B, FitInfo] = lasso(X, y, 'CV', cv);
                    save shit.mat
                    idx = FitInfo.Index1SE;
                    coef = B(:, idx);
                    coef0 = FitInfo.Intercept(idx);
                    pred = data(s).betas{c} * coef + coef0;
                    mse(s) = immse(y, pred(~data(s).exclude));

                otherwise
                    assert(false);
            end

            pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
            decRU = [decRU; pred];
        end
        exclude = logical(exclude);

        tbl = data2table(data, 0, 0); % include all trials; we exclude bad runs and timeouts manually
        tbl = [tbl table(decRU)];

        
        % glm with both RU and decRU
        results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
        [w, names, stats] = fixedEffects(results_both{c});
        ps(c,:) = stats.pValue';
        results_both{c}
        stats.pValue
        w
        names

        % glm with RU only
        % do model comparison
        results_RU{c} = fitglme(tbl,formula_RU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
        comp{c} = compare(results_RU{c}, results_both{c}); % order is important -- see docs
        comp{c}
        p_comp(c,:) = comp{c}.pValue(2);
        BIC(c,:) = comp{c}.BIC';

        % glm with decRU only
        % do second model comparison
        results_decRU{c} = fitglme(tbl,formula_decRU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
        comp2{c} = compare(results_decRU{c}, results_both{c}); % order is important -- see docs
        comp2{c}
        p_comp2(c,:) = comp2{c}.pValue(2);
        BIC2(c,:) = comp2{c}.BIC';

        % correlate RMSE with behavioral weights across subjects
        % => see if better decodeability is associated with more reliance on regressor in decision
        %
        load results_glme_fig3_nozscore.mat;
        w = getEffects(results_VTURU, false);
        [r, p] = corr(abs(w(:,2)), mse');
        disp('mse to w');
        r
        p
        p_ax(c,:) = p;
        r_ax(c,:) = r;
    end

    save(filename, '-v7.3');

    p_uncorr = ps(:,4);
    p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
    BIC_RU = BIC(:,1);
    BIC_both = BIC(:,2);
    BIC_decRU = BIC2(:,1);
    disp(method);
    table(masknames', p_uncorr, p_corr, BIC_RU, BIC_both, p_comp, BIC_decRU, p_comp2, p_ax, r_ax)

end



function pred = ridgepred(X, y, Xtest, Lambda)
    coef = ridge(y, X, Lambda, 0);
    Xtest = [ones(size(Xtest, 1), 1), Xtest]; % include intercept term
    pred = Xtest * coef;
end
