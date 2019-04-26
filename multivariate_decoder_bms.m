% multivariate decoder analysis 
% streamlined version of multivariate_decoder, with BMS
%
% see if activation in ROI predicts choices better than regressor from model
%
function multivariate_decoder_bms(roi_glmodel, roi_contrast, regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth)

printcode;

rmpath('/n/sw/helmod/apps/centos7/Core/spm/12.7487-fasrc01/external/fieldtrip/external/stats/'); % for binopdf on cluster

assert(standardize ~= 1, 'Don''t z-score! It makes the w''s meaningless, also it''s incorrect.');

nTRs = 242;
[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

data = load_data;

null_iters = 100;

if ~exist('do_orth', 'var')
    do_orth = false;
end
if ~exist('standardize', 'var')
    standardize = false;
end
if ~exist('mixed_effects', 'var')
    mixed_effects = false;
end
if ~exist('clusterFWEcorrect', 'var')
    clusterFWEcorrect = true;
end
if ~exist('extent', 'var')
    extent = [];
end
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end
if ~exist('intercept', 'var')
    intercept = false; 
end
if ~exist('zscore_across_voxels', 'var')
    zscore_across_voxels = false; % whether to z-score betas across voxels
end
if ~exist('predict_abs', 'var')
    predict_abs = false; % whether to predict |RU| instead of RU, |V| instead of V, etc; flip_sign = true by default (see univariate_decoder)
end
if ~exist('use_smooth', 'var')
    use_smooth = false; % whether to use smooth activations
end


filename = sprintf('multivariate_decoder_bms_roiglm%d_%s_%s_orth=%d_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d_intercept=%d_method=%s_getnull=%d_zav=%d_pa=%d_us=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth);
disp(filename);

% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num);


% define behavioral / hybrid GLM formulas
[formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept);
formula_both
formula_orig



beta_series_glm = 23;

if use_smooth
    EXPT = exploration_expt();
else
    EXPT = exploration_expt_nosmooth();
end


% extract behavioral regressors
% e.g. make them match the rows in the data table (NaNs for missing runs)
% or flip RU based on sign, etc
%
V_all = [];
DV_all = [];
exclude = [];
for s = 1:length(data)
    [V, RU, TU, VTU, DV] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    data(s).V = V;
    data(s).DV = DV;
    V_all = [V_all; V];
    DV_all = [DV_all; DV];

    runs = find(goodRuns{s}); % only those runs were included in the GLMs
    data(s).bad_runs = ~ismember(data(s).run, runs); % ... those runs were NOT included in the GLMs
    data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials

    exclude = [exclude; data(s).exclude];
end
exclude = logical(exclude);


save(filename, '-v7.3');

best_of = 3; % get best model (BIC-wise) out of how many

% original behavioral glm   

tbl = data2table(data, standardize, 0); % include all trials; we exclude bad runs and timeouts manually

successes = 0;
results_orig = [];
for attempt = 1:100
    try
        res = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', 'random');
        res 

        if isempty(results_orig) || results_orig.LogLikelihood < res.LogLikelihood
            results_orig = res;
        end
        successes = successes + 1;

        if successes == best_of
            break;
        end
    catch e
        fprintf('             failed fitting "%s" on attempt %d...\n', formula_orig, attempt);
        disp(e)
    end
end
assert(attempt < 100, 'failed too many times...');

[BICs, logliks] = get_subj_bics(results_orig, tbl, exclude);
disp('Original behavioral GLM');
results_orig

LMEs = [-0.5 * BICs];


% fit behavioral GLM with activations
%
ps = [];
null_ps = [];
for c = 1:numel(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);

    fprintf('Region: %s\n', region{c});

    dec = [];
    mse = [];
    null_p = [];
    for s = 1:length(data)
        % get beta series
        B = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', masks{c});
        assert(size(B,1) > 1);
        assert(size(B,2) > 1);

        % exclude nan voxels
        which_nan = any(isnan(B), 1); % exclude nan voxels (ignoring bad runs and timeouts; we exclude those in the GLMs)
        B(:, which_nan) = [];

        % init betas with # trials = # of rows in behavioral data
        X = nan(length(data(s).run), size(B,2));

        % not all runs were used in GLM => only set betas for the trials we will use
        X(~data(s).bad_runs,:) = B;

        if zscore_across_voxels
            % optionally get rid of mean BOLD signal
            X = zscore(X,[],2);
        end

        switch regressor
            case 'RU'
                if predict_abs
                    y = abs(data(s).RU);
                else
                    y = data(s).RU;
                end
            case 'TU'
                y = data(s).TU;
            case 'V'
                if predict_abs
                    y = abs(data(s).V);
                else
                    y = data(s).V;
                end
            case 'DV'
                if predict_abs
                    y = abs(data(s).DV);
                else
                    y = data(s).DV;
                end
            otherwise
                assert(false);
        end

        % predict, excluding bad trials (timeouts & bad runs) 
        % for CV, one run per fold
        [pred, mse(s)] = multilinear_fit(X, y, X, method, data(s).run, data(s).exclude);
        data(s).mse{c} = mse(s);

        % pad up predictions with nan's for bad trials 
        tmp = pred;
        pred = nan(length(data(s).run), 1);
        pred(~data(s).exclude) = tmp;

        if predict_abs
            % account for fact that we're predicing e.g. |RU| and not RU
            % same as flip_sign in univariate_decoder
            switch regressor
                case 'RU'
                    pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0);
                case 'V'
                    pred = pred .* (data(s).V >= 0) + (-pred) .* (data(s).V < 0);
                case 'DV'
                    pred = pred .* (data(s).DV >= 0) + (-pred) .* (data(s).DV < 0);
            end
        end

        dec = [dec; pred];
        data(s).dec{c} = pred;

        assert(sum(isnan(data(s).dec{c}(~data(s).exclude))) == 0, 'got NaN predictions');

        % optionally generate null distribution
        if get_null
            null_mse = [];
            for i = 1:null_iters
                y = y(randperm(length(y)));
                [~, m] = multilinear_fit(X, y, X, method, data(s).run, data(s).exclude);
                null_mse = [null_mse, m];
            end
            data(s).null_mse{c} = null_mse;

            % calculate p-value based on null distribution
            null_mse = [null_mse mse(s)];
            null_mse = sort(null_mse);
            idx = find(null_mse == mse(s));
            p = idx(1) / length(null_mse);
            fprintf('                    subj %d null mse p = %.4f\n', s, p);
            data(s).null_p{c} = p;
            null_p(s) = p;
        end
    end

    if get_null
        null_ps = [null_ps; null_p];
    end

    tbl_dec = augment_table_with_decoded_regressor(tbl, regressor, dec, standardize, exclude, V_all);

    %
    % fitglme sometimes gets NaN log likelihood and fails, especially for random effects 
    % => need to try a few times with random starts
    %

    % augmented glm with decoded regressor 
    successes = 0;
    results_both{c} = [];
    BIC(c,:) = [NaN NaN]; % in case it doesn't work
    p_comp(c,:) = NaN;
    for attempt = 1:100
        try
            res = fitglme(tbl_dec,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', 'random');
            [w, names, stats] = fixedEffects(res);
            res
            stats.pValue
            w

            assert(res.LogLikelihood >= results_orig.LogLikelihood, 'Loglik of augmented model is no better than original model');

            save wtf.mat
            if isempty(results_both{c}) || results_both{c}.LogLikelihood < res.LogLikelihood
                results_both{c} = res;

                ps(c,:) = stats.pValue';
                comp{c} = compare(results_orig, results_both{c}); % order is important -- see docs
                comp{c}
                p_comp(c,:) = comp{c}.pValue(2);
                BIC(c,:) = comp{c}.BIC';

                [~, reg_names] = fixedEffects(results_both{c});
            end
            successes = successes + 1;

            if successes == best_of
                break;
            end
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula_orig, attempt);
            disp(e)
        end
    end

    if successes > 0
        BICs = get_subj_bics(results_both{c}, tbl_dec, exclude);
        LMEs = [LMEs, -0.5 * BICs];
    end


    % sanity check -- activations should correlate with regressor
    switch regressor
        case 'RU'
            RU = table2array(tbl(:,'RU'));
            [r,p] = corr(RU(~exclude), dec(~exclude));
        case 'TU'
            TU = table2array(tbl(:,'TU'));
            [r,p] = corr(TU(~exclude), dec(~exclude));
        case 'V'
            V = table2array(tbl(:,'V'));
            [r,p] = corr(V(~exclude), dec(~exclude));
        case 'DV'
            [r,p] = corr(DV_all(~exclude), dec(~exclude));
    end


    pears_rs(c,:) = r;
    pears_ps(c,:) = p;

    % correlate MSE with behavioral weights across subjects
    % => see if better decodeability is associated with more reliance on regressor in decision
    %
    if standardize == 1
        assert(false, 'NO ZSCORING!');
        load results_glme_fig3.mat; % TODO does not exist
    elseif standardize == 2
        load results_glme_fig3_norm.mat;
    else
        load results_glme_fig3_nozscore.mat;
    end

    w = getEffects(results_VTURU, false);
    switch regressor
        case 'RU'
            [r, p] = corr(abs(w(:,2)), mse');
        case 'TU'
            [r, p] = corr(abs(w(:,3)), mse');
        case 'V'
            [r, p] = corr(abs(w(:,1)), mse');
        case 'DV'
            r = NaN;
            p = NaN;
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

[alpha,exp_r,xp,pxp,bor] = bms(LMEs);

fprintf('BOR = %.6f\n', bor);
fprintf('PXP of original GLM = %.6f\n', pxp(1));
pxp = pxp(2:end);
pxp = pxp';

reg_idx = find(contains(reg_names.Name, 'dec'));

p_uncorr = ps(:,reg_idx);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
p_comp_corr = 1 - (1 - p_comp) .^ numel(p_comp);


if get_null
    frac_s = mean(null_ps < 0.05, 2); % fraction of participants whose null distribution p-value is < 0.05, i.e. we can significantly decode regressor
    table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, p_comp_corr, pxp, p_ax, r_ax, frac_s)
else
    table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, p_comp_corr, pxp, p_ax, r_ax)
end

