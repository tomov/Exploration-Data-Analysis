% univariate decoder analysis 
% see if activation in ROI predicts choices better than regressor from model
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function univariate_decoder(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept)

printcode;

assert(standardize ~= 1, 'Don''t z-score! It makes the w''s meaningless, also it''s incorrect.');

nTRs = 242;
EXPT = exploration_expt();
[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

data = load_data;

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

filename = sprintf('univariate_decoder_roiglm%d_%s_glm%d_%s_orth=%d_lambda=%f_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d_intercept=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept);
disp(filename);

% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num);


% find closest TR to each trial onset (adjusted for HRF f'n)
for s = 1:length(data)
    runs = find(goodRuns{s}); % only those runs were included in the GLMs
    data(s).bad_runs = ~ismember(data(s).run, runs); % ... those runs were NOT included in the GLMs

    % trial onset idx's from bad runs are NaNs
    [~,session] = find(data(s).run == runs); % automatically excludes bad runs
    data(s).trial_onset_act_idx = nan(length(data(s).trial_onset), 1);
    data(s).trial_onset_act_idx(~data(s).bad_runs) = get_activations_idx(EXPT, data(s).trial_onset(~data(s).bad_runs), session, nTRs);
end


% define behavioral / hybrid GLM formulas
[formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept);
formula_both
formula_orig


% decode regressor 
for c = 1:length(masks)
    for s = 1:length(data)
        dec = ccnl_decode_regressor(EXPT, glmodel, regressor, masks{c}, lambda, s);
        data(s).all_act{c} = dec{1};
    end
end

% massage decoded regressors
% e.g. make them match the rows in the data table (NaNs for missing runs)
% or flip RU based on sign, etc
%
V_all = [];
DV_all = [];
for s = 1:length(data)
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
    load(fullfile(modeldir,'SPM.mat'));

    data(s).act = nan(length(data(s).run), length(masks));
    [V, RU, TU, VTU, DV] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    data(s).V = V;
    data(s).DV = DV;
    V_all = [V_all; V(~data(s).timeout)];
    DV_all = [DV_all; DV(~data(s).timeout)];

    for c = 1:length(masks)
        % pick trial_onset activations only
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations (excluding bad runs, which were excluded in the GLM)
        act{c} = data(s).all_act{c}(which_act,:); % only consider 1 activation for each trial

        % average across voxels in ROI
        % notice not all runs were used in the GLMs
        data(s).act(~data(s).bad_runs,c) = mean(act{c}, 2);

        % adjust for fact that the regressor was |RU|
        %if strcmp(regressor, 'RU')
        %    data(s).act(:,c) = data(s).act(:,c) .* (RU >= 0) + (-data(s).act(:,c)) .* (RU < 0);
        %end

        %% adjust for fact that the regressor was |V|
        %if strcmp(regressor, 'V')
        %    data(s).act(:,c) = data(s).act(:,c) .* (V >= 0) + (-data(s).act(:,c)) .* (V < 0);
        %end

        %% adjust for fact that the regressor was |DV|
        %if strcmp(regressor, 'DV')
        %    data(s).act(:,c) = data(s).act(:,c) .* (DV >= 0) + (-data(s).act(:,c)) .* (DV < 0);
        %end
    end
end

save(filename, '-v7.3');


% fit behavioral GLM with activations
%
ps = [];
for c = 1:numel(masks)
    act = [];
    mse = [];
    bad_runs = logical([]);
    for s = 1:length(data)
        act = [act; data(s).act(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        bad_runs = [bad_runs; data(s).bad_runs(~data(s).timeout)]; % bad runs are also out (their activations are NaNs)

        which = ~data(s).bad_runs & ~data(s).timeout;
        switch regressor % TODO act is still whitened & filtered => MSE might be wrong
            case 'RU'
                mse(s) = immse(data(s).RU(which), data(s).act(which, c));
            case 'TU'
                mse(s) = immse(data(s).TU(which), data(s).act(which, c));
            case 'V'
                mse(s) = immse(data(s).V(which), data(s).act(which, c));
            case 'DV'
                mse(s) = immse(data(s).DV(which), data(s).act(which, c));
        end
    end
    assert(all(isnan(act(bad_runs))));
    assert(all(~isnan(act(~bad_runs))));

    tbl = data2table(data,standardize,1); % exclude timeouts for fitting

    tbl = augment_table_with_decoded_regressor(tbl, regressor, act, standardize, bad_runs, V_all);

    %
    % fitglme sometimes gets NaN log likelihood and fails, especially for random effects 
    % => need to try a few times with random starts
    %

    % augmented glm with decoded regressor 
    for attempt = 1:100
        try
            results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',bad_runs, 'StartMethod', 'random', 'verbose', 2);
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
            results_orig{c} = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace','CovariancePattern','diagonal','EBMethod','TrustRegion2D', 'Exclude',bad_runs, 'StartMethod', 'random', 'verbose', 2);
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
            [r,p] = corr(RU(~bad_runs), act(~bad_runs));
        case 'TU'
            TU = table2array(tbl(:,'TU'));
            [r,p] = corr(TU(~bad_runs), act(~bad_runs));
        case 'V'
            V = table2array(tbl(:,'V'));
            [r,p] = corr(V(~bad_runs), act(~bad_runs));
        case 'DV'
            [r,p] = corr(DV_all(~bad_runs), act(~bad_runs));
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


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
p_comp_corr = 1 - (1 - p_comp) .^ numel(p_comp);
table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, p_comp_corr, p_ax, r_ax)

