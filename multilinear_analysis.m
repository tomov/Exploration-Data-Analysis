% multilinear regression analysis 
% try to decode |RU| or TU from multivariate ROI activity and see if it predicts
% choices better than behavioral model alone
% merge of residuals_analysis.m and badre_2012_multilinear_analysis.m TODO dedupe?
%

function multilinear_analysis(glmodel, regressor, contrast, load_first_half)

printcode;

if ~exist('load_first_half', 'var')
    load_first_half = false;
end

filename = ['multilinear_analysis_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '.mat'];
filename


if load_first_half
    % optionally load pre-computed residuals (b/c the 2nd half doesn't work on the stupid cluster...)
    load(filename);
else

    % group-level settings
    p = 0.001;
    alpha = 0.05;
    Dis = 20;
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
    direct = '+';

    EXPT = exploration_expt();

    data = load_data;

    switch regressor
        case 'RU'
            formula_both = 'C ~ -1 + V + RU + VTU + decRU';
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + decRU + VTU';
        case 'TU'
            formula_both = 'C ~ -1 + V + RU + VTU + VdecTU';
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + RU + VdecTU';
        otherwise
            assert(false);
    end

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

    save(filename, '-v7.3');
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

for c = 1:numel(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);
    disp(region{c});

    dec = [];
    exclude = [];
    rmse = [];
    for s = 1:length(data)
        exclude = [exclude; data(s).exclude];
        X = data(s).betas{c};
        y = data(s).y;
        mdl = fitlm(X, y, 'exclude', data(s).exclude, 'Intercept', true);
        pred = predict(mdl, X);
        if glmodel == 21 && strcmp(regressor, 'RU')
            pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
        end
        dec = [dec; pred];
        rmse(s) = mdl.RMSE;
    end
    exclude = logical(exclude);

    tbl = data2table(data, 0, 0); % include all trials; we exclude bad runs and timeouts manually

    switch regressor
        case 'RU'
            decRU = dec;
            tbl = [tbl table(decRU)];
        case 'TU'
            VdecTU = V_all ./ dec;
            tbl = [tbl table(VdecTU)];
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

    % correlate RMSE with behavioral weights across subjects
    % => see if better decodeability is associated with more reliance on regressor in decision
    %
    load results_glme_fig3_nozscore.mat;
    w = getEffects(results_VTURU, false);
    switch regressor
        case 'RU'
            [r, p] = corr(abs(w(:,2)), rmse');
        case 'TU'
            [r, p] = corr(abs(w(:,3)), rmse');
        otherwise
            assert(false);
    end
    disp('rmse to w');
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

