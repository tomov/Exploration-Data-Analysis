% multivariate decoder analysis 
% merge of univariate_decoder and multilinear_analysis
%
% see if activation in ROI predicts choices better than regressor from model
%
function multivariate_decoder(roi_glmodel, roi_contrast, regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, method, get_null)

printcode;

assert(standardize ~= 1, 'Don''t z-score! It makes the w''s meaningless, also it''s incorrect.');

nTRs = 242;
EXPT = exploration_expt();
[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

data = load_data;

betas_from_mat = false;
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


filename = sprintf('multivariate_decoder_roiglm%d_%s_%s_orth=%d_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d_intercept=%d_method=%s_getnull=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, method, get_null);
disp(filename);

% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num);


% define behavioral / hybrid GLM formulas
[formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept);
formula_both
formula_orig


% --------- extract betas from .mat files -----------

if betas_from_mat
    % extract betas (GLM 57, saved as .mat files)
    roi = extract_roi_betas(masks, 'trial_onset');

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
end

% ---------------------------------------------


if ~betas_from_mat
    beta_series_glm = 57;
    EXPT_nosmooth = exploration_expt_nosmooth();
end

% extract & massage betas
% e.g. make them match the rows in the data table (NaNs for missing runs)
% or flip RU based on sign, etc
%
V_all = [];
DV_all = [];
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

    % --------- alternatively, extract betas from disk TODO sanity check -----------
    %
    if ~betas_from_mat
        for c = 1:length(masks)
            % get beta series
            B = get_beta_series(EXPT_nosmooth, beta_series_glm, s, 'trial_onset', masks{c});

            % exclude nan voxels
            which_nan = any(isnan(B), 1); % exclude nan voxels (ignoring bad runs and timeouts; we exclude those in the GLMs)
            B(:, which_nan) = [];

            % init betas with # trials = # of rows in behavioral data
            data(s).betas{c} = nan(length(data(s).run), size(B,2));

            % not all runs were used in GLM => only set betas for the trials we will use
            data(s).betas{c}(~data(s).bad_runs,:) = B;
        end
    end
end


save(filename, '-v7.3');


% fit behavioral GLM with activations
%
ps = [];
for c = 1:numel(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);
    disp(region{c});

    dec = [];
    exclude = [];
    mse = [];
    for s = 1:length(data)
        X = data(s).betas{c};

        switch regressor
            case 'RU'
                y = data(s).RU;
            case 'TU'
                y = data(s).TU;
            case 'V'
                y = data(s).V;
            case 'DV'
                y = data(s).DV;
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

        % sometimes predictions are nan's => exclude those too
        exclude = [exclude; data(s).exclude | isnan(pred)];

        %if strcmp(regressor, 'RU')
        %    pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
        %end NB: we don't do this any more b/c we directly predict RU
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

    tbl = data2table(data, standardize, 0); % include all trials; we exclude bad runs and timeouts manually

    tbl = augment_table_with_decoded_regressor(tbl, regressor, dec, standardize, exclude, V_all);

    %
    % fitglme sometimes gets NaN log likelihood and fails, especially for random effects 
    % => need to try a few times with random starts
    %

    % augmented glm with decoded regressor 
    for attempt = 1:100
        try
            results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', 'random', 'verbose', 2);
            [w, names, stats] = fixedEffects(results_both{c});
            ps(c,:) = stats.pValue';
            results_both{c}
            stats.pValue
            w
            break
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula_both, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times');


    % original glm  
    % do model comparison
    for attempt = 1:100
        try
            results_orig{c} = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', 'random', 'verbose', 2);
            comp{c} = compare(results_orig{c}, results_both{c}); % order is important -- see docs
            comp{c}
            p_comp(c,:) = comp{c}.pValue(2);
            BIC(c,:) = comp{c}.BIC';
            break
        catch e
            fprintf('             failed fitting "%s" on attempt %d...\n', formula_orig, attempt);
            disp(e)
        end
    end
    assert(attempt < 100, 'failed too many times');


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


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
p_comp_corr = 1 - (1 - p_comp) .^ numel(p_comp);
table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, p_comp_corr, p_ax, r_ax)

