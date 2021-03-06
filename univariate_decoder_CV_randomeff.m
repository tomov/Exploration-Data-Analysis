% univariate decoder analysis from univarite_decoder_refactored.m, but with CV like neurosynth_CV.m (TODO dedupe)
% see if activation in ROI predicts choices better than regressor from model
%
% pick based lambda based on independent GLMs for each other subject (rather than mixed effects) => doesn't really work, probs not enough power in each subject
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m
% TODO dedupe w/ multivariate_decoder_bms and univariate_decoder

function univariate_decoder_CV(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, get_null)

printcode;

rmpath('/n/sw/helmod/apps/centos7/Core/spm/12.7487-fasrc01/external/fieldtrip/external/stats/'); % for binopdf on cluster

assert(standardize ~= 1, 'Don''t z-score! It makes the w''s meaningless, also it''s incorrect.');

nTRs = 242;
EXPT = exploration_expt();
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
if ~exist('flip_sign', 'var')
    flip_sign = false; % flip the sign of the decoded |RU|, |V|, |DV|, etc based on the sign of RU
end
if ~exist('get_null', 'var')
    get_null = false;  % generate null distribution
end

Lambda = logspace(-10,10,21);

GLM_has_timeouts = true; % does glmodel include timeout trials?
assert(GLM_has_timeouts, 'sorry this is a hardcoded assumption');

filename = sprintf('univariate_decoder_CV_randomeff_roiglm%d_%s_glm%d_%s_orth=%d_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d_intercept=%d_flip=%d_gn=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, get_null);
disp(filename);


% define behavioral / hybrid GLM formulas
[formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept);
formula_both
formula_orig


% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num);


% extract behavioral regressors & stuff
%
V_all = [];
DV_all = [];
exclude = [];
for s = 1:length(data)

    % get behavioral regressors
    [V, RU, TU, VTU, DV] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    data(s).V = V;
    data(s).DV = DV;
    V_all = [V_all; V];
    DV_all = [DV_all; DV];

    % figure out which trials to exclude from hybrid GLMs (timeouts & bad runs with too much motion)
    runs = find(goodRuns{s}); % only those runs were included in the GLMs
    data(s).bad_runs = ~ismember(data(s).run, runs); % ... those runs were NOT included in the GLMs
    data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials

    exclude = [exclude; data(s).exclude];

    % initialize empty decoded array
    data(s).act = nan(length(data(s).run), length(Lambda), length(masks));

    % figure out which TRs toughly correspond to trial_onsets (with HRF offset & stuff)
    %
    % trial onset idx's from bad runs are NaNs
    % note that timeouts are included here b/c we do have them in the fMRI GLM (not the behavioral one though)
    [~,session] = find(data(s).run == runs); % automatically excludes bad runs
    data(s).trial_onset_act_idx = nan(length(data(s).trial_onset), 1);
    data(s).trial_onset_act_idx(~data(s).bad_runs) = get_activations_idx(EXPT, data(s).trial_onset(~data(s).bad_runs), session, nTRs);
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

BICs = get_subj_bics(results_orig, tbl, exclude);
disp('Original behavioral GLM');
results_orig

LMEs = [-0.5 * BICs];

formula_both_fixed = remove_random_effects(formula_both);

% fit behavioral GLM with activations
%
null_ps = [];
ps = [];
for c = 1:numel(masks)

    disp(masks{c});

    logliks{c} = []; % analogous to mses in neurosynth_CV

    % first pass -- compute (random effects) hybrid GLM loglik for different lambdas
    %
    null_p = [];
    for s = 1:length(data)
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
        load(fullfile(modeldir,'SPM.mat'));

        % decode regressor
        dec = ccnl_decode_regressor(EXPT, glmodel, ['x', regressor, '^'], masks{c}, Lambda, s);
        dec = mean(dec{1}, 2); % average across voxels
        dec = squeeze(dec); % so its [TRs, lambdas]

        % pick trial_onset activations only
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations (excluding bad runs, which were excluded in the GLM)

        % subset activations & trials
        % on the left, we only assign to runs that were used in the GLM
        % on the right, we only take decoded regressor from trial onsets (+ 5 sec HRF offset)
        data(s).act(~data(s).bad_runs,:,c) = dec(which_act,:);

        % optionally flip sign
        data(s).act(:,:,c) = adjust_sign(data(s), data(s).act(:,:,c), regressor, flip_sign);

        % for each lambda, fit augmented GLM (for that subject only) and get loglik
        %
        tbl_s = data2table(data(s), standardize, 0);

        for l = 1:length(Lambda)

            tbl_s_dec = augment_table_with_decoded_regressor(tbl_s, regressor, data(s).act(:,l,c), standardize, data(s).exclude, data(s).V);

            % assume won't fail b/c no random effects
            res_s = fitglme(tbl_s_dec, formula_both_fixed, 'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',data(s).exclude, 'StartMethod', 'default'); 

            data(s).loglik(l,c) = res_s.LogLikelihood;
            logliks{c}(s,l) = res_s.LogLikelihood;
        end

        % optionally generate null distribution
        if get_null
            assert(false, 'unsupported'); 
        end
    end

    if get_null
        null_ps = [null_ps; null_p];
    end


    % second pass -- decode using best lambda from other subjects
    %

    act = [];
    mse = [];
    for s = 1:length(data)
        % find lambda that gives best total loglik for other subjects
        total_loglik = sum(logliks{c}([1:s-1 s+1:end],:), 1); % analogous to avg_mse in neurosynth_CV
        %[~, idx] = max(total_loglik);
        idx = 11; % TODO sanity

        Lambda_idxs{c}(s) = idx;

        fprintf('                        subj %d -- lambda(%d) = %f\n', s, idx, Lambda(idx));

        % use that lambda for decoding
        act = [act; data(s).act(:, idx, c)];

        mse(s) = calc_mse(data(s), data(s).act(:, idx, c), regressor, flip_sign);
    end
    assert(all(~isnan(act(~exclude))));


    tbl_dec = augment_table_with_decoded_regressor(tbl, regressor, act, standardize, exclude, V_all);

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
            if attempt == successes + 1
                StartMethod = 'default'; % prefer default start method, unless it's failing on us
            else
                StartMethod = 'random'; % if it's failing, try random start method
            end
            StartMethod

            res = fitglme(tbl_dec,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',exclude, 'StartMethod', StartMethod);
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
            if flip_sign
                RU = tbl.RU;
            else
                RU = abs(tbl.RU);
            end
            [r,p] = corr(RU(~exclude), act(~exclude));
        case 'TU'
            TU = tbl.TU;
            [r,p] = corr(TU(~exclude), act(~exclude));
        case 'V'
            if flip_sign
                V = tbl.V; 
            else
                V = abs(tbl.V);
            end
            [r,p] = corr(V(~exclude), act(~exclude));
        case 'DV'
            if flip_sign
                DV = tbl.DV;
            else
                DV = abs(tbl.DV);
            end
            [r,p] = corr(DV_all(~exclude), act(~exclude));
    end


    pears_rs(c,:) = r;
    pears_ps(c,:) = p;

    % correlate MSE with behavioral weights across subjects
    % => see if better decodeability is associated with more reliance on regressor in decision
    %
    if standardize == 1
        assert(false, 'NO ZSCORING!');
        load results_glme_fig3_TrustRegion2D.mat; % TODO does not exist
    elseif standardize == 2
        load results_glme_fig3_norm_TrustRegion2D.mat;
    else
        load results_glme_fig3_nozscore_TrustRegion2D.mat;
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
if size(pxp,1) == 1
    pxp = pxp';
end

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

end


% pass data(s)
%
function mse = calc_mse(data, dec, regressor, flip_sign)
    which = ~data.bad_runs & ~data.timeout;
    switch regressor % TODO act is still whitened & filtered => MSE might be wrong
        case 'RU'
            if flip_sign
                mse = immse(data.RU(which), dec(which));
            else
                mse = immse(abs(data.RU(which)), dec(which));
            end
        case 'TU'
            mse = immse(data.TU(which), dec(which));
        case 'V'
            if flip_sign
                mse = immse(data.V(which), dec(which));
            else
                mse = immse(abs(data.V(which)), dec(which));
            end
        case 'DV'
            if flip_sign
                mse = immse(data.DV(which), dec(which));
            else
                mse = immse(abs(data.DV(which)), dec(which));
            end
    end
end


% pass data(s)
%
function dec = adjust_sign(data, dec, regressor, flip_sign)
    % flip the sign of the decoded regressor
    % while this may seem like double dipping / cheating, we are including RU, V, etc in the GLM, so this sign flip does not add new information
    if flip_sign
        switch regressor
            case 'RU'
                % adjust for fact that the regressor was |RU|
                dec = dec .* (data.RU >= 0) + (-dec) .* (data.RU < 0);
            case 'V'

                % adjust for fact that the regressor was |V|
                dec = dec .* (data.V >= 0) + (-dec) .* (data.V < 0);

            case 'DV'
                % adjust for fact that the regressor was |DV|
                dec = dec .* (data.DV >= 0) + (-dec) .* (data.DV < 0);
        end
    end
end

% remove random effects component from fitglme formula
%
function formula = remove_random_effects(formula)

    st = strfind(formula, '+ (');
    if isempty(st) % no random effects component
        return; 
    end
    assert(length(st) == 1);

    en = strfind(formula, ')');
    assert(length(en) == 1);

    formula = erase(formula, formula(st:en));
end
