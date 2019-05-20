% multivariate decoder analysis, whole-brain: neurosynth parcellation (merge of multivariate_decoder_bms and is_rsa)
% streamlined version of multivariate_decoder, with BMS
%
% see if activation in ROI predicts choices better than regressor from model
%
function neurosynth_CV(regressor, do_orth, standardize, mixed_effects, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth, lateralized, parcel_idxs)

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
if ~exist('lateralized', 'var')
    lateralized = false;
end

% get ROIs
group_mask_filename = fullfile('masks', 'mask.nii');

if lateralized
    parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2_lateralized.nii');
else
    parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2.nii');
end

[~, Vparcel, parcellation_vol] = ccnl_load_mask(parcellation_file);
parcellation_vol = round(parcellation_vol);

if ~exist('parcel_idxs', 'var') || isempty(parcel_idxs)
    parcel_idxs = unique(parcellation_vol(:));
else
    assert(all(ismember(parcel_idxs, unique(parcellation_vol(:)))));
end


filename = sprintf('neurosynth_CV_%s_orth=%d_standardize=%d_mixed=%d_intercept=%d_method=%s,%s_getnull=%d_zav=%d_pa=%d_us=%d_l=%d_pi=%d.mat', regressor, do_orth, standardize, mixed_effects, intercept, method{1}, method{2}, get_null, zscore_across_voxels, predict_abs, use_smooth, lateralized, length(parcel_idxs));
disp(filename);


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


% extract behavioral regressors & stuff
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

%
% original behavioral glm   
%

tbl = data2table(data, standardize, 0); % include all trials; we exclude bad runs and timeouts manually

successes = 0;
results_orig = [];
for attempt = 1:100
    try
        if attempt == successes + 1
            StartMethod = 'default'; % prefer default start method, unless it's failing on us
        else
            StartMethod = 'random'; % if it's failing, try random start method
        end
        StartMethod

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

Lambda = logspace(-10,20,31);

% fit behavioral GLM with activations
%
ps = [];
c = 0;
region = [];
null_ps = [];
for i = 1:length(parcel_idxs)

    parcel_idx = parcel_idxs(i);
    if parcel_idx == 0
        continue;
    end

    i
    region = [region; parcel_idx];
    fprintf('parcel = %d\n', parcel_idx);

    % ROI
    mask = parcellation_vol == parcel_idx;

    % normalize mask
    group_vol = spm_vol(group_mask_filename);
    group_mask = spm_read_vols(group_vol);

    [x, y, z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
    cor = mni2cor(cor2mni([x y z], Vparcel.mat), group_vol.mat); % voxel coords in AAL2 space --> voxel coords in MNI space --> voxel coords in our space
    ind = sub2ind(size(group_mask), cor(:,1), cor(:,2), cor(:,3)); % voxel coords in our space --> voxel indices
    
    % Reconstruct mask in our space
    %
    Vmask = group_vol;
    Vmask.fname = 'tmp.nii'; % CRUCIAL! o/w overwrite mask.nii
    mask = zeros(size(group_mask));
    mask(ind) = 1; % voxel indices --> binary mask
    
    % Only include voxels that are part of the subject group-level mask
    % i.e. that have non-NaN betas for all subjects
    %
    mask = mask & group_mask;

    if sum(mask(:)) == 0
        disp('skipping parcel -- empty mask');
        continue; % some masks are already lateralized
    end

    % first pass -- compute MSE for different lambdas
    %

    % ok we're sure we're using this mask
    c = c + 1;

    fprintf('      first pass\n');

    mses{c} = [];

    for s = 1:length(data)
        fprintf('              subj %d -- getting betas\n', s);

        % get beta series
        tic
        B = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', mask);
        toc
        assert(size(B,1) > 1);
        assert(size(B,2) == sum(mask(:)));

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

        % only take mses
        fprintf('              subj %d -- computing fit for all lambdas\n', s);
        tic
        [~, ~, mses{c}(s,:)] = multilinear_fit(X, y, X, method{1}, data(s).run, data(s).exclude, Lambda);
        toc

        % hopefully speed things up for the second pass
        all_B{s} = B;
        all_X{s} = X;
        all_y{s} = y;
    end


    % second pass -- decode using best lambda from other subjects
    %

    dec = [];
    mse = [];
    null_p = [];
    for s = 1:length(data)
        % reuse from 1st pass 
        B = all_B{s};
        X = all_X{s};
        y = all_y{s};

        % find lambda that gives best average MSE for other subjects
        avg_mse = mean(mses{c}([1:s-1 s+1:end],:), 1);
        [~, idx] = min(avg_mse);
        Lambda_idxs{c}(s) = idx;

        fprintf('                        subj %d -- lambda(%d) = %f\n', s, idx, Lambda(idx));

        % predict, excluding bad trials (timeouts & bad runs) 
        % for CV, one run per fold
        tic
        [pred, mse(s)] = multilinear_fit(X, y, X, method{2}, data(s).run, data(s).exclude, Lambda(idx));
        toc
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
            runs = unique(data(s).run);
            for iter = 1:null_iters
                % shuffle runs !!! non-exchangeability at single trial level, even at block level
                % b/c they're dependent
                runs_perm = runs(randperm(length(runs)));
                for r = 1:length(runs_perm)
                    y_perm(data(s).run == runs(r),:) = y(data(s).run == runs_perm(r),:);
                end
                [~, m] = multilinear_fit(X, y_perm, X, method{2}, data(s).run, data(s).exclude, Lambda(idx));
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

    decs{c} = dec;

    tbl_dec = augment_table_with_decoded_regressor(tbl, regressor, dec, standardize, exclude, V_all);

    fprintf('                          fitting GLMs for parcel %d\n', parcel_idx);
    tic

    %
    % fitglme sometimes gets NaN log likelihood and fails, especially for random effects 
    % => need to try a few times with random starts
    %

    % augmented glm with decoded regressor 
    successes = 0;
    results_both = [];
    BIC(c,:) = [NaN NaN]; % in case it doesn't work
    p_comp(c,:) = NaN;
    for attempt = 1:100
        try
            if attempt == successes + 1
                StartMethod = 'default'; % prefer default start method, unless it's failing on us
            else
                StartMethod = 'random'; % if it's failing, try random start method
            end
            StartMethod
            
            res = fitglme(tbl_dec,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', 'random');
            [w, names, stats] = fixedEffects(res);
            res
            stats.pValue
            w

            assert(res.LogLikelihood >= results_orig.LogLikelihood, 'Loglik of augmented model is no better than original model');

            if isempty(results_both) || results_both.LogLikelihood < res.LogLikelihood
                results_both = res;

                ps(c,:) = stats.pValue';
                comp{c} = compare(results_orig, results_both); % order is important -- see docs
                comp{c}
                p_comp(c,:) = comp{c}.pValue(2);
                BIC(c,:) = comp{c}.BIC';

                [~, reg_names] = fixedEffects(results_both);
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

    toc

    if successes > 0
        BICs = get_subj_bics(results_both, tbl_dec, exclude);
        LMEs = [LMEs, -0.5 * BICs];
    end


    % sanity check -- activations should correlate with regressor
    % IMPORTANT -- correlate within each run (fold) separately
    % then do t-test on correlation coefficients
    rs = [];
    ps = [];
    for s = 1:length(data)
        for r = 1:max(data(s).run)
            which = ~exclude & tbl.S == s & tbl.run == r;
            if sum(which) == 0
                continue; % bad run
            end
            switch regressor
                case 'RU'
                    RU = table2array(tbl(:,'RU'));
                    [r,p] = corr(RU(which), dec(which));
                case 'TU'
                    TU = table2array(tbl(:,'TU'));
                    [r,p] = corr(TU(which), dec(which));
                case 'V'
                    V = table2array(tbl(:,'V'));
                    [r,p] = corr(V(which), dec(which));
                case 'DV'
                    [r,p] = corr(DV_all(which), dec(which));
            end
            rs = [rs r];
            ps = [ps p];
        end
    end

    pears_rs_all{c} = rs;
    pears_ps_all{c} = ps;

    pears_rs(c,:) = mean(rs); % take mean r before Fisher z transform

    rs = atanh(rs);
    [h, p, ci, stat] = ttest(rs);

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


    fprintf('                  saving to %s\n', filename);

    tic
    save(filename, '-v7.3'); % save after every parcel, just in case
    toc
end



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
pears_ps_corr = 1 - (1 - pears_ps) .^ numel(pears_ps);

if get_null
    frac_s = mean(null_ps < 0.05, 2); % fraction of participants whose null distribution p-value is < 0.05, i.e. we can significantly decode regressor
    table(region, p_uncorr, p_corr, pears_rs, pears_ps, pears_ps_corr, BIC_orig, BIC_both, p_comp, p_comp_corr, pxp, p_ax, r_ax, frac_s)
else
    table(region, p_uncorr, p_corr, pears_rs, pears_ps, pears_ps_corr, BIC_orig, BIC_both, p_comp, p_comp_corr, pxp, p_ax, r_ax)
end
