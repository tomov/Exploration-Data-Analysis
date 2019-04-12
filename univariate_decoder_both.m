% univariate decoder analysis for both RU & TU
% see if activation in ROI predicts choices better than regressor from model
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function univariate_decoder_both(glmodel, RU_roi_idx, TU_roi_idx, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent)

printcode;

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

filename = sprintf('univariate_decoder_both_glm%d_RUroi=%d_TUroi=%d_orth=%d_lambda=%f_standardize=%d_mixed=%d_corr=%d_extent=%d.mat', glmodel, RU_roi_idx, TU_roi_idx, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent);
disp(filename);

% get ROIs
[masks_RU, region_RU] = get_masks(glmodel, 'RU', clusterFWEcorrect, extent);
[masks_TU, region_TU] = get_masks(glmodel, 'TU', clusterFWEcorrect, extent);
masks{1} = masks_RU{RU_roi_idx};
masks{2} = masks_TU{TU_roi_idx};
region{1,:} = region_RU{RU_roi_idx};
region{2,:} = region_TU{TU_roi_idx};

regressor = {'RU', 'TU'};


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
[formula_both, formula_orig] = get_formula('both', do_orth, mixed_effects);


% decode regressor 
for c = 1:length(masks)
    for s = 1:length(data)
        dec = ccnl_decode_regressor(EXPT, glmodel, regressor{c}, masks{c}, lambda, s);
        data(s).all_act{c} = dec{1};
    end
end

% massage decoded regressors
% e.g. make them match the rows in the data table (NaNs for missing runs)
% or flip RU based on sign, etc
%
V_all = [];
for s = 1:length(data)
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
    load(fullfile(modeldir,'SPM.mat'));

    data(s).act = nan(length(data(s).run), length(masks));
    [V, RU, TU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    V_all = [V_all; V(~data(s).timeout)];

    for c = 1:length(masks)
        % not all runs were used in the GLMs
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations (excluding bad runs, which were excluded in the GLM)
        act{c} = data(s).all_act{c}(which_act,:); % only consider 1 activation for each trial

        % average across voxels in ROI
        data(s).act(~data(s).bad_runs,c) = mean(act{c}, 2);

        % adjust for fact that the regressor was |RU|
        if strcmp(regressor{c}, 'RU')
            data(s).act(:,c) = data(s).act(:,c) .* (RU >= 0) + (-data(s).act(:,c)) .* (RU < 0);
        end
    end
end

save(filename, '-v7.3');

% massage activations some more
%
for c = 1:numel(masks)
    act{c} = [];
    mse{c} = [];
    bad_runs = logical([]);
    for s = 1:length(data)
        act{c} = [act{c}; data(s).act(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        bad_runs = [bad_runs; data(s).bad_runs(~data(s).timeout)]; % bad runs are also out (their activations are NaNs)

        which = ~data(s).bad_runs & ~data(s).timeout;
        switch regressor{c} % TODO act is still whitened & filtered => MSE might be wrong
            case 'RU'
                mse{c}(s) = immse(data(s).RU(which), data(s).act(which, c));
            case 'TU'
                mse{c}(s) = immse(data(s).TU(which), data(s).act(which, c));
        end
    end
    assert(all(isnan(act{c}(bad_runs))));
    assert(all(~isnan(act{c}(~bad_runs))));
end

% get model-based regressors
%
tbl = data2table(data,standardize,1); % exclude timeouts for fitting

% add decoded regressors to tbl
%
for c = 1:numel(masks)
    switch regressor{c}
        case 'RU'
            decRU = act{c};
            if standardize == 1
                decRU(~bad_runs) = zscore(decRU(~bad_runs));
            elseif standardize == 2
                decRU(~bad_runs) = decRU(~bad_runs) / norm(decRU(~bad_runs));
            end
            tbl = [tbl table(decRU)];

            % orthogonalized version
            tmp = spm_orth([tbl.RU(~bad_runs), decRU(~bad_runs)]);
            decRU_orth = decRU;
            decRU_orth(~bad_runs) = tmp(:,2);
            if standardize == 1
                decRU_orth(~bad_runs) = zscore(decRU_orth(~bad_runs));
            elseif standardize == 2
                decRU_orth(~bad_runs) = decRU_orth(~bad_runs) / norm(decRU_orth(~bad_runs));
            end
            tbl = [tbl table(decRU_orth)];

        case 'TU'
            VdecTU = V_all ./ act{c};
            if standardize == 1
                VdecTU(~bad_runs) = zscore(VdecTU(~bad_runs));
            elseif standardize == 2
                VdecTU(~bad_runs) = VdecTU(~bad_runs) / norm(VdecTU(~bad_runs));
            end
            tbl = [tbl table(VdecTU)];

            % orthogonalized version
            tmp = spm_orth([tbl.VTU(~bad_runs), VdecTU(~bad_runs)]);
            VdecTU_orth = VdecTU;
            VdecTU_orth(~bad_runs) = tmp(:,2); 
            if standardize == 1
                VdecTU_orth(~bad_runs) = zscore(VdecTU_orth(~bad_runs));
            elseif standardize == 2
                VdecTU_orth(~bad_runs) = VdecTU_orth(~bad_runs) / norm(VdecTU_orth(~bad_runs));
            end
            tbl = [tbl table(VdecTU_orth)];
        otherwise
            assert(false);
    end

    % sanity check -- activations should correlate with regressor
    switch regressor{c}
        case 'RU'
            RU = table2array(tbl(:,'RU'));
            [r,p] = corr(RU(~bad_runs), act{c}(~bad_runs));
        case 'TU'
            TU = table2array(tbl(:,'TU'));
            [r,p] = corr(TU(~bad_runs), act{c}(~bad_runs));
    end

    pears_rs(c,:) = r;
    pears_ps(c,:) = p;
end


% fit behavioral GLM with activations
%
ps = [];

% glm with both RU/TU and actRU/actTU
results_both = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',bad_runs);
[w, names, stats] = fixedEffects(results_both);
results_both
stats.pValue
w

% glm with RU/TU only
% do model comparison
results_orig = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',bad_runs);
comp = compare(results_orig, results_both); % order is important -- see docs
comp
%p_comp(c,:) = comp.pValue(2);
%BIC(c,:) = comp.BIC';



results_all = results_both;

%load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat', 'results_both');
load('univariate_decoder_roiglm36_RU_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat', 'results_both');
results_RU = results_both{RU_roi_idx};

%load('univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat', 'results_both');
load('univariate_decoder_roiglm36_TU_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat', 'results_both');
results_TU = results_both{TU_roi_idx};

comp_RU = compare(results_RU, results_all);
comp_TU = compare(results_TU, results_all);

comp_RU

comp_TU

save(filename, '-v7.3');

%{

p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
BIC_dec = BIC2(:,1);
table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, BIC_dec, p_comp2, p_ax, r_ax)
%}
