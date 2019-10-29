% copy of univariate_decoder_refactored but for correlating residuals (i.e. (Vhat - V)^2 with TU^2)
% TODO dedupe with univariate_decoder_refactored and all stuff in there

function univariate_decoder_residuals(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, do_CV, get_null, sphere)

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
if ~exist('lambda', 'var')
    lambda = 1;
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
if ~exist('do_CV', 'var')
    do_CV = false;  % cross-validate betas, i.e. beta for each run = avg beta from all other runs
end
if ~exist('get_null', 'var')
    get_null = false;  % generate null distribution
end
if ~exist('sphere', 'var')
    sphere = 10; % 10 mm sphere by default
end

GLM_has_timeouts = true; % does glmodel include timeout trials?
assert(GLM_has_timeouts, 'sorry this is a hardcoded assumption');


filename = sprintf('univariate_decoder_residuals_roiglm%d_%s_glm%d_%s_orth=%d_lambda=%f_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d_intercept=%d_flip=%d_doCV=%d_gn=%d_s=%.1f.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, do_CV, get_null, sphere);
disp(filename);


% define behavioral / hybrid GLM formulas
[formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept);
formula_both
formula_orig


% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num, sphere);
masks'


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
    data(s).VTU = VTU;
    data(s).DV = DV;
    V_all = [V_all; V];
    DV_all = [DV_all; DV];

    % figure out which trials to exclude from hybrid GLMs (timeouts & bad runs with too much motion)
    runs = find(goodRuns{s}); % only those runs were included in the GLMs
    data(s).bad_runs = ~ismember(data(s).run, runs); % ... those runs were NOT included in the GLMs
    data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials

    exclude = [exclude; data(s).exclude];

    % initialize empty decoded array
    data(s).act = nan(length(data(s).run), length(masks));

    % figure out which TRs toughly correspond to trial_onsets (with HRF offset & stuff)
    %
    % trial onset idx's from bad runs are NaNs
    % note that timeouts are included here b/c we do have them in the fMRI GLM (not the behavioral one though)
    [~,session] = find(data(s).run == runs); % automatically excludes bad runs
    data(s).trial_onset_act_idx = nan(length(data(s).trial_onset), 1);
    data(s).trial_onset_act_idx(~data(s).bad_runs) = get_activations_idx(EXPT, data(s).trial_onset(~data(s).bad_runs), session, nTRs);
end
exclude = logical(exclude);





% fit behavioral GLM with activations
%
null_ps = [];
ps = [];
for c = 1:numel(masks)

    disp(masks{c});

    % TODO first pass -- compute (random effects) hybrid GLM loglik for different lambdas
    %
    null_p = [];
    for s = 1:length(data)
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
        load(fullfile(modeldir,'SPM.mat'));

        % decode regressor
        dec = ccnl_decode_regressor(EXPT, glmodel, ['x', regressor, '^'], masks{c}, lambda, s);
        dec = mean(dec{1}, 2); % average across voxels

        % pick trial_onset activations only
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations (excluding bad runs, which were excluded in the GLM)

        % subset activations & trials
        % on the left, we only assign to runs that were used in the GLM
        % on the right, we only take decoded regressor from trial onsets (+ 5 sec HRF offset)
        data(s).act(~data(s).bad_runs,c) = dec(which_act,:);

        % optionally flip sign
        data(s).act(:,c) = adjust_sign(data(s), data(s).act(:,c), regressor, flip_sign);

        % calculate residual = (y - yhat)
        assert(flip_sign); % assume sign is correct

        switch regressor
            case 'RU'
                y = data(s).RU;

            case 'V'
                y = data(s).V;

            case 'VTU'
                y = data(s).VTU;

            case 'DV'
                y = data(s).DV;

            otherwise
                assert(false);
        end

        y = y(~data(s).bad_runs); % model V
        yhat = data(s).act(~data(s).bad_runs,c); % Vhat, decoded from the brain
        v = data(s).TU(~data(s).bad_runs);

        [r,p] = corr((y - yhat).^2, v.^2); % TODO glm not corr TODO sqrt?

        data(s).corr_r(c) = r;
        data(s).corr_p(c) = p;
        corr_rs(c,s) = r;

    end


    % TODO second pass -- decode using best lambda from other subjects
    %

    act = [];
    mse = [];
    for s = 1:length(data)
        act = [act; data(s).act(:, c)];

        mse(s) = calc_mse(data(s), data(s).act(:,c), regressor, flip_sign);
    end
    assert(all(~isnan(act(~exclude))));


    tbl_dec = augment_table_with_decoded_regressor(tbl, regressor, act, standardize, exclude, V_all);

end


corr_zrs = atanh(corr_rs);

for c = 1:numel(masks)
    [h,p,ci,stat] = ttest(corr_zrs(c,:));

    ts(c,:) = stat.tstat;
    ps(c,:) = p;
    dfs(c,:) = stat.df;
end


save(filename, '-v7.3');

table(region, ts, dfs, ps)


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
        case 'VTU'
            if flip_sign
                mse = immse(data.VTU(which), dec(which));
            else
                mse = immse(abs(data.VTU(which)), dec(which));
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

            case 'VTU'
                % adjust for fact that the regressor was |VTU|
                dec = dec .* (data.VTU >= 0) + (-dec) .* (data.VTU < 0);

            case 'DV'
                % adjust for fact that the regressor was |DV|
                dec = dec .* (data.DV >= 0) + (-dec) .* (data.DV < 0);
        end
    end
end
