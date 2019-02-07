% univariate decoder analysis 
% see if activation in ROI predicts choices better than regressor from model
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function univariate_decoder(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num)

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
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end

filename = sprintf('univariate_decoder_roiglm%d_%s_glm%d_%s_orth=%d_lambda=%f_standardize=%d_mixed=%d_corr=%d_extent=%d_Num=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num);
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
[formula_both, formula_orig, formula_dec] = get_formula(regressor, do_orth, mixed_effects);


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
for s = 1:length(data)
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
    load(fullfile(modeldir,'SPM.mat'));

    data(s).act = nan(length(data(s).run), length(masks));
    [V, RU, TU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    data(s).V = V;
    V_all = [V_all; V(~data(s).timeout)];

    for c = 1:length(masks)
        % not all runs were used in the GLMs
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations (excluding bad runs, which were excluded in the GLM)
        act{c} = data(s).all_act{c}(which_act,:); % only consider 1 activation for each trial

        % average across voxels in ROI
        data(s).act(~data(s).bad_runs,c) = mean(act{c}, 2);

        % adjust for fact that the regressor was |RU|
        if strcmp(regressor, 'RU')
            data(s).act(:,c) = data(s).act(:,c) .* (RU >= 0) + (-data(s).act(:,c)) .* (RU < 0);
        end

        % adjust for fact that the regressor was |V|
        if strcmp(regressor, 'V')
            data(s).act(:,c) = data(s).act(:,c) .* (V >= 0) + (-data(s).act(:,c)) .* (V < 0);
        end
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
        end
    end
    assert(all(isnan(act(bad_runs))));
    assert(all(~isnan(act(~bad_runs))));

    tbl = data2table(data,standardize,1); % exclude timeouts for fitting

    % TODO dedupe with univariate_decoder_both and three
    switch regressor
        case 'RU'
            decRU = act;
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
            VdecTU = V_all ./ act;
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

        case 'V'
            decV = act;
            if standardize == 1
                decV(~bad_runs) = zscore(decV(~bad_runs));
            elseif standardize == 2
                decV(~bad_runs) = decV(~bad_runs) / norm(decV(~bad_runs));
            end
            tbl = [tbl table(decV)];

            % orthogonalized version
            tmp = spm_orth([tbl.V(~bad_runs), decV(~bad_runs)]);
            decV_orth = decV;
            decV_orth(~bad_runs) = tmp(:,2);
            if standardize == 1
                decV_orth(~bad_runs) = zscore(decV_orth(~bad_runs));
            elseif standardize == 2
                decV_orth(~bad_runs) = decV_orth(~bad_runs) / norm(decV_orth(~bad_runs));
            end
            tbl = [tbl table(decV_orth)];


        otherwise
            assert(false);
    end

    % glm with both RU and actRU
    results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',bad_runs);
    [w, names, stats] = fixedEffects(results_both{c});
    ps(c,:) = stats.pValue';
    results_both{c}
    stats.pValue
    w

    % glm with RU only
    % do model comparison
    results_orig{c} = fitglme(tbl,formula_orig,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',bad_runs);
    comp{c} = compare(results_orig{c}, results_both{c}); % order is important -- see docs
    comp{c}
    p_comp(c,:) = comp{c}.pValue(2);
    BIC(c,:) = comp{c}.BIC';

    % glm with actRU only
    % do second model comparison
    results_dec{c} = fitglme(tbl,formula_dec,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',bad_runs);
    comp2{c} = compare(results_dec{c}, results_both{c}); % order is important -- see docs
    comp2{c}
    p_comp2(c,:) = comp2{c}.pValue(2);
    BIC2(c,:) = comp2{c}.BIC';


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
    end


    pears_rs(c,:) = r;
    pears_ps(c,:) = p;

    % correlate MSE with behavioral weights across subjects
    % => see if better decodeability is associated with more reliance on regressor in decision
    %
    if standardize == 1
        load results_glme_fig3.mat;
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
table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, BIC_dec, p_comp2, p_ax, r_ax)

