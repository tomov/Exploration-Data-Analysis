% univariate decoder analysis 
% see if activation in ROI predicts choices better than regressor from model
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function univariate_decoder(glmodel, regressor, contrast, normalize, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent)

printcode;

EXPT = exploration_expt();

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

filename = sprintf('univariate_decoder_glm%d_%s_%s_norm=%d_orth=%d_lambda=%f_standardize=%d_mixed=%d_corr=%d_extent=%d.mat', glmodel, regressor, replace(contrast, ' ', '_'), normalize, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent);
disp(filename);

% get ROI masks
switch contrast
    case 'badre'
        % clusters = masks from paper
        masks = badre_2012_create_masks(false);
        %masks = masks(1); % TODO use all masks

        for c = 1:length(masks)
            mask = masks{c};
            [~, masknames{c}, ~] = fileparts(mask);
            region{c,:} = masknames{c};
        end


    case 'dlpfc'
        % clusters = masks from paper
        masks = dlpfc_2012_create_masks(false);

        for c = 1:length(masks)
            mask = masks{c};
            [~, masknames{c}, ~] = fileparts(mask);
            region{c,:} = masknames{c};
        end


    case 'tommy'
        % clusters = masks from paper
        masks = tommy_2017_create_masks(false);

        for c = 1:length(masks)
            mask = masks{c};
            [~, masknames{c}, ~] = fileparts(mask);
            region{c,:} = masknames{c};
        end

    otherwise
        % group-level settings
        p = 0.001;
        alpha = 0.05;
        Dis = 20;
        Num = 1; % # peak voxels per cluster; default in bspmview is 3
        direct = '+';

        [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

        r = 10 / 1.5; % 10 mm radius

        % create spherical masks around peak voxel of each cluster (intersected with cluster)
        %
        for c = 1:length(region)
            masks{c} = sprintf('sphere_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));
            cmask = CI == CI(cor(c,1), cor(c,2), cor(c,3));
            ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c}, cmask);
        end

end

% find peak of HRF
hrf = spm_hrf(0.001);
[~,hrf_offset] = max(hrf);
hrf_offset = hrf_offset / 1000;

nTRs = 242;
TR = EXPT.TR;
trs = TR/2 : TR : nTRs * TR;

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

% find closest TR to each trial onset (adjusted for HRF f'n)
for s = 1:length(data)
    act_idx = [];
    runs = find(goodRuns{s});
    data(s).bad_runs = ~ismember(data(s).run, runs); % exclude bad runs
    for i = 1:length(data(s).trial_onset)
        [~, idx] = min(abs(trs - (data(s).trial_onset(i) + hrf_offset)));
        if data(s).bad_runs(i)
            act_idx = [act_idx; NaN];
        else
            r = find(data(s).run(i) == runs); % scan session idx in GLM 
            act_idx = [act_idx; idx + nTRs * (r - 1)];
        end
    end
    data(s).trial_onset_act_idx = act_idx;
end


% define behavioral / hybrid GLM formulas
switch regressor
    case 'RU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + (-1 + V + RU + VTU + decRU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + (-1 + V + RU + VTU + decRU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
            formula_dec = 'C ~ -1 + V + decRU + VTU + (-1 + V + decRU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + decRU + VTU';
        end

    case 'TU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth + (-1 + V + RU + VTU + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU + (-1 + V + RU + VTU + VdecTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
            formula_dec = 'C ~ -1 + V + RU + VdecTU + (-1 + V + RU + VdecTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + RU + VdecTU';
        end

    otherwise
        assert(false);
end


% get betas to (optionally) normalize activations in each run
for c = 1:length(masks)
    mask = masks{c};
    m = load_mask(mask);
    cnt = sum(m(:));

    for s = 1:length(data)
        runs = find(goodRuns{s});
        data(s).b{c} = nan(length(data(s).run), cnt);

        if normalize == 0
            % do nothing
        elseif normalize == 1
            % act_RU = (act - b0) / b_RU
            % i.e. assume other b's are insignificant
            %
            for run = 1:max(data(s).run)
                r = find(run == runs); % scan session idx in GLM
                if ~isempty(r)
                    % get beta for regressor
                    reg = ['Sn(', num2str(r), ') trial_onsetx', regressor];
                    fprintf('  c = %s, s = %d, run = %d, r = %d, reg = %s\n', mask, s, run, r, reg);
                    b = ccnl_get_beta(EXPT, glmodel, reg, mask, s);
                    data(s).b{c}(data(s).run == run, :) = repmat(b, sum(data(s).run == run), 1);

                    % get beta0
                    reg = ['Sn(', num2str(r), ') constant'];
                    fprintf('  c = %s, s = %d, run = %d, r = %d, reg = %s\n', mask, s, run, r, reg);
                    b0 = ccnl_get_beta(EXPT, glmodel, reg, mask, s);
                    data(s).b0{c}(data(s).run == run, :) = repmat(b0, sum(data(s).run == run), 1);
                end
            end
        elseif normalize == 2
            % act_RU = (act - X_\RU * b_\RU) ./ b_RU
            % i.e. take other regressors into accoutn
            %
            % TODO dedupe with ccnl_get_beta and ccnl_get_activations / ccnl_get_residuals
            % also improve those based on this
            %
            modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
            load(fullfile(modeldir,'SPM.mat'));
            names = SPM.xX.name';
            cdir = pwd;
            cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
            B = spm_data_read(SPM.Vbeta, find(m));
            cd(cdir);
            X = SPM.xX.X;

            % separate RU betas and regressors from the rest
            which_reg = contains(names, regressor);
            B_noreg = B(~which_reg, :);
            B_reg = B(which_reg, :);
            B_reg = repelem(B_reg, nTRs, 1); % we're need one for each TR b/c we're doing element-wise divison by b_RU
            X_noreg = X(:, ~which_reg);
            X_reg = X(:, which_reg);

            fprintf('  c = %s, s = %d\n', mask, s);

            data(s).B_noreg{c} = B_noreg;
            data(s).B_reg{c} = B_reg;
            data(s).X_noreg{c} = X_noreg;
            data(s).X_reg{c} = X_reg;

        elseif normalize == 3 || normalize == 4
            % act_RU = (act - X_\RU * b_\RU) ./ b_RU
            % i.e. take other regressors into accoutn
            % same as 2 BUT using whitened / filtered X and Y (!) like SPM
            %
            modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
            load(fullfile(modeldir,'SPM.mat'));
            names = SPM.xX.name';
            cdir = pwd;
            cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
            B = spm_data_read(SPM.Vbeta, find(m));
            cd(cdir);
            X = SPM.xX.xKXs.X;

            % separate RU betas and regressors from the rest
            which_reg = contains(names, regressor);
            B_noreg = B(~which_reg, :);
            B_reg = B(which_reg, :);
            B_reg = repelem(B_reg, nTRs, 1); % we're need one for each TR b/c we're doing element-wise divison by b_RU
            X_noreg = X(:, ~which_reg);
            X_reg = X(:, which_reg);

            fprintf('  c = %s, s = %d\n', mask, s);

            data(s).B_noreg{c} = B_noreg;
            data(s).B_reg{c} = B_reg;
            data(s).X_noreg{c} = X_noreg;
            data(s).X_reg{c} = X_reg;
        else
            assert(false);

        end
    end
end

% extract activations for each cluster
%
V_all = [];
for s = 1:length(data)
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
    load(fullfile(modeldir,'SPM.mat'));

    clear act;
    for c = 1:length(masks)
        mask = masks{c};
        [~, masknames{c}, ~] = fileparts(mask);

        act{c} = ccnl_get_activations(EXPT, glmodel, mask, s);
        data(s).all_act{c} = act{c};

    end

    data(s).act = nan(length(data(s).run), length(masks));
    [V, RU, TU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    data(s).RU = RU;
    data(s).TU = TU;
    V_all = [V_all; V(~data(s).timeout)];

    for c = 1:length(masks)
        if normalize == 2
            act{c} = (act{c} - data(s).X_noreg{c} * data(s).B_noreg{c}) ./ data(s).B_reg{c};
        elseif normalize == 3
            act{c} = spm_filter(SPM.xX.K,SPM.xX.W*act{c});
            act{c} = (act{c} - data(s).X_noreg{c} * data(s).B_noreg{c}) ./ data(s).B_reg{c};
        elseif normalize == 4
            % ridge regression -- regulalize by lambda
            % x_RU = (activation - sum of x_i * b_i, for i != RU) * b_RU / (b_RU^2 + lambda)
            % strictly speaking we should call it decRU instead of act but whatevs
            %
            act{c} = spm_filter(SPM.xX.K,SPM.xX.W*act{c});
            act{c} = (act{c} - data(s).X_noreg{c} * data(s).B_noreg{c}) .* data(s).B_reg{c} ./ (data(s).B_reg{c}.^2 + lambda);
        end

        % not all runs were used in the GLMs
        which_act = data(s).trial_onset_act_idx(~data(s).bad_runs); % trial onset activations
        act{c} = act{c}(which_act,:); % only consider 1 activation for each trial

        if normalize == 1
            act{c} = (act{c} - data(s).b0{c}(~data(s).bad_runs)) ./ data(s).b{c}(~data(s).bad_runs);
        end

        data(s).act(~data(s).bad_runs,c) = mean(act{c}, 2);

        % adjust for fact that the regressor was |RU|
        if glmodel == 21 && strcmp(regressor, 'RU')
            data(s).act(:,c) = data(s).act(:,c) .* (RU >= 0) + (-data(s).act(:,c)) .* (RU < 0);
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
        end
    end
    assert(all(isnan(act(bad_runs))));
    assert(all(~isnan(act(~bad_runs))));

    tbl = data2table(data,standardize,1); % exclude timeouts for fitting

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

