% multilinear regression analysis 
% try to decode |RU| or TU from multivariate ROI activity and see if it predicts
% choices better than behavioral model alone
% merge of residuals_analysis.m and badre_2012_multilinear_analysis.m TODO dedupe?
%

function multilinear_analysis(glmodel, regressor, contrast, method, get_null, do_orth, load_first_half)

printcode;

null_iters = 100;

if ~exist('get_null', 'var')
    get_null = false;
end

if ~exist('do_orth', 'var')
    do_orth = true;
end

if ~exist('load_first_half', 'var')
    load_first_half = false;
end

filename = ['multilinear_analysis_', regressor, '_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '_', method, '_orth=', num2str(do_orth), '.mat'];
filename

if load_first_half
    % optionally load pre-computed betas to speed things up
    load(filename);
else

    EXPT = exploration_expt();

    data = load_data;

    % special case a priori ROIs; otherwise expect a glm & contrast
    switch glmodel
        case 'badre'
            % clusters = masks from paper
            masks = badre_2012_create_masks(false);
            %masks = masks(1); % TODO all masks
            region = masks';

            % extract trial_onset (raw, unsmoothed) betas
            roi = extract_roi_betas(masks, 'trial_onset');

        case 'tommy'
            % clusters = masks from paper
            masks = tommy_2017_create_masks(false);
            region = masks';

            % extract trial_onset (raw, unsmoothed) betas
            roi = extract_roi_betas(masks, 'trial_onset');

        otherwise

            % group-level settings
            p = 0.001;
            alpha = 0.05;
            Dis = 20;
            Num = 1; % # peak voxels per cluster; default in bspmview is 3
            direct = '+';

            [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);

            r = 10 / 1.5; % 10 mm radius

            % create spherical masks around peak voxel of each cluster (intersected with cluster)
            %
            for c = 1:length(region)
                masks{c} = sprintf('sphere_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));
                cmask = CI == CI(cor(c,1), cor(c,2), cor(c,3));
                ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c}, cmask);
            end

            % extract trial_onset (raw, unsmoothed) betas
            %
            roi = extract_roi_betas(masks, 'trial_onset');
    end

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


    % define behavioral / hybrid GLM formulas
    switch regressor
        case 'RU'
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + decRU + VTU';

        case 'TU'
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + RU + VdecTU';

        otherwise
            assert(false);
    end

    % extract regressors
    %
    V_all = [];
    for s = 1:length(data)
        which_all = logical(ones(length(data(s).run), 1));
        switch regressor
            case 'RU'
                [~, absRU] =  get_latents(data, s, which_all, 'abs');
                [~, RU] = get_latents(data, s, which_all, 'left');
                data(s).y = absRU;
                data(s).RU = RU; % for sign-correction
            case 'TU'
                [V, ~, TU] = get_latents(data, s, which_all, 'left');
                data(s).y = TU;
                V_all = [V_all; V];
            otherwise
                assert(false);
        end
    end

    save(filename, '-v7.3');
end

rng default; % for reproducibility

% run analysis
%
for c = 1:numel(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);
    disp(region{c});

    dec = [];
    exclude = [];
    mse = [];
    for s = 1:length(data)
        exclude = [exclude; data(s).exclude];
        X = data(s).betas{c};
        y = data(s).y;

        % remove bad data points
        X = X(~data(s).exclude, :);
        y = y(~data(s).exclude);

        % predict using full data set; we ignore bad trials later 
        % also for CV, one run per fold
        [pred, mse(s)] = multilinear_fit(X, y, data(s).betas{c}, method, data(s).run(~data(s).exclude));
        data(s).mse{c} = mse(s);

        if strcmp(method, 'ridge_CV_CV')
            tmp = pred;
            pred = nan(length(data(s).run), 1);
            pred(~data(s).exclude) = tmp;
        end

        if strcmp(regressor, 'RU')
            pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
        end
        dec = [dec; pred];

        % optionally generate null distribution
        if get_null
            null_mse = [];
            for i = 1:null_iters
                y = y(randperm(length(y)));
                [~, m] = multilinear_fit(X, y, data(s).betas{c}, method, data(s).run(~data(s).exclude));
                null_mse = [null_mse, m];
            end
            data(s).null_mse{c} = null_mse;

            % calculate p-value based on null distribution
            null_mse = [null_mse mse(s)];
            null_mse = sort(null_mse);
            idx = find(null_mse == mse(s));
            p = idx / length(null_mse);
            fprintf('                    subj %d null mse p = %.4f\n', s, p);
            data(s).null_p{c} = p;
        end
    end
    exclude = logical(exclude);

    tbl = data2table(data, 0, 0); % include all trials; we exclude bad runs and timeouts manually

    switch regressor
        case 'RU'
            decRU = dec;
            tbl = [tbl table(decRU)];
            % orthogonalized version
            tmp = spm_orth([tbl.RU(~exclude), decRU(~exclude)]);
            decRU_orth = decRU;
            decRU_orth(~exclude) = tmp(:,2);
            tbl = [tbl table(decRU_orth)];
        case 'TU'
            VdecTU = V_all ./ dec;
            tbl = [tbl table(VdecTU)];
            % orthogonalized version
            tmp = spm_orth([tbl.VTU(~exclude), VdecTU(~exclude)]);
            VdecTU_orth = VdecTU;
            VdecTU_orth(~exclude) = tmp(:,2);
            tbl = [tbl table(VdecTU_orth)];
        otherwise
            assert(false);
    end

    % glm with both the model and decoded regressor
    results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results_both{c});
    ps(c,:) = stats.pValue';
    results_both{c}
    stats.pValue
    w
    names

    % original glm (model only)
    % do model comparison
    results_orig{c} = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp{c} = compare(results_orig{c}, results_both{c}); % order is important -- see docs
    comp{c}
    p_comp(c,:) = comp{c}.pValue(2);
    BIC(c,:) = comp{c}.BIC';

    % glm with decoded regressor only
    % do second model comparison
    results_dec{c} = fitglme(tbl,formula_dec,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp2{c} = compare(results_dec{c}, results_both{c}); % order is important -- see docs
    comp2{c}
    p_comp2(c,:) = comp2{c}.pValue(2);
    BIC2(c,:) = comp2{c}.BIC';

    % correlate MSE with behavioral weights across subjects
    % => see if better decodeability is associated with more reliance on regressor in decision
    %
    load results_glme_fig3_nozscore.mat;
    w = getEffects(results_VTURU, false);
    switch regressor
        case 'RU'
            [r, p] = corr(abs(w(:,2)), mse');
        case 'TU'
            [r, p] = corr(abs(w(:,3)), mse');
        otherwise
            assert(false);
    end
    disp('mse to w');
    r
    p
    p_ax(c,:) = p;
    r_ax(c,:) = r;
end

save(filename, '-v7.3');

p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
BIC_dec = BIC2(:,1);
table(region, p_uncorr, p_corr, BIC_orig, BIC_both, p_comp, BIC_dec, p_comp2, p_ax, r_ax)

