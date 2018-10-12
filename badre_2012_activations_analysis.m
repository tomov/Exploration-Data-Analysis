% activations analysis for RU for Badre 2012 RLPFC ROI
% see if average activation in ROI predicts choices better than RU from model
%
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function badre_2012_activations_analysis(glmodel, normalize)

printcode;

EXPT = exploration_expt();

if ~exist('glmodel', 'var')
    glmodel = 21;
end
if ~exist('normalize', 'var')
    normalize = 1; % divide each activation by the corresponding beta
end

data = load_data;

formula_both = 'C ~ -1 + V + RU + VTU + actRU';
formula_RU = 'C ~ -1 + V + RU + VTU';
formula_actRU = 'C ~ -1 + V + actRU + VTU';

filename = ['badre_2012_activations_analysis_glm', num2str(glmodel), '_normalize', num2str(normalize), '.mat'];
disp(filename);

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
    data(s).exclude = ~ismember(data(s).run, runs); % exclude bad runs
    for i = 1:length(data(s).trial_onset)
        [~, idx] = min(abs(trs - (data(s).trial_onset(i) + hrf_offset)));
        if data(s).exclude(i)
            act_idx = [act_idx; NaN];
        else
            r = find(data(s).run(i) == runs); % scan session idx in GLM 
            act_idx = [act_idx; idx + nTRs * (r - 1)];
        end
    end
    data(s).trial_onset_act_idx = act_idx;
end


% clusters = masks from paper
masks = badre_2012_create_masks(false);
masks = masks(1); % TODO use all masks

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
                    % get beta for RU
                    reg = ['Sn(', num2str(r), ') trial_onsetxRU'];
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
            which_RU = contains(names, 'RU');
            B_noRU = B(~which_RU, :);
            B_RU = B(which_RU, :);
            B_RU = repelem(B_RU, nTRs, 1); % we're need one for each TR b/c we're doing element-wise divison by b_RU
            X_noRU = X(:, ~which_RU);
            X_RU = X(:, which_RU);

            fprintf('  c = %s, s = %d\n', mask, s);

            data(s).B_noRU{c} = B_noRU;
            data(s).B_RU{c} = B_RU;
            data(s).X_noRU{c} = X_noRU;
            data(s).X_RU{c} = X_RU;

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
            which_RU = contains(names, 'RU');
            B_noRU = B(~which_RU, :);
            B_RU = B(which_RU, :);
            B_RU = repelem(B_RU, nTRs, 1); % we're need one for each TR b/c we're doing element-wise divison by b_RU
            X_noRU = X(:, ~which_RU);
            X_RU = X(:, which_RU);

            fprintf('  c = %s, s = %d\n', mask, s);

            data(s).B_noRU{c} = B_noRU;
            data(s).B_RU{c} = B_RU;
            data(s).X_noRU{c} = X_noRU;
            data(s).X_RU{c} = X_RU;
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
    [V, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    V_all = [V_all; V(~data(s).timeout)];

    for c = 1:length(masks)
        if normalize == 2
            act{c} = (act{c} - data(s).X_noRU{c} * data(s).B_noRU{c}) ./ data(s).B_RU{c};
        elseif normalize == 3
            act{c} = spm_filter(SPM.xX.K,SPM.xX.W*act{c});
            act{c} = (act{c} - data(s).X_noRU{c} * data(s).B_noRU{c}) ./ data(s).B_RU{c};
        elseif normalize == 4
            % ridge regression -- regulalize by lambda
            % x_RU = (activation - sum of x_i * b_i, for i != RU) * b_RU / (b_RU^2 + lambda)
            % strictly speaking we should call it decRU instead of act but whatevs
            %
            lambda = 1; % TODO determine in principled way; hard to do b/c GLM is already run => nothing to fit / cross-validate...
            act{c} = spm_filter(SPM.xX.K,SPM.xX.W*act{c});
            act{c} = (act{c} - data(s).X_noRU{c} * data(s).B_noRU{c}) .* data(s).B_RU{c} ./ (data(s).B_RU{c}.^2 + lambda);
        end

        % not all runs were used in the GLMs
        which_act = data(s).trial_onset_act_idx(~data(s).exclude); % trial onset activations
        act{c} = act{c}(which_act,:); % only consider 1 activation for each trial

        if normalize == 1
            act{c} = (act{c} - data(s).b0{c}(~data(s).exclude)) ./ data(s).b{c}(~data(s).exclude);
        end

        data(s).act(~data(s).exclude,c) = mean(act{c}, 2);

        % adjust for fact that the regressor was |RU|
        if glmodel == 21
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
    exclude = logical([]);
    for s = 1:length(data)
        act = [act; data(s).act(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        exclude = [exclude; data(s).exclude(~data(s).timeout)]; % bad runs are also out (their activations are NaNs)
    end
    assert(all(isnan(act(exclude))));
    assert(all(~isnan(act(~exclude))));

    tbl = data2table(data,0,1); % exclude timeouts for fitting
    actRU = act;
    tbl = [tbl table(actRU)];

    % glm with both RU and actRU
    results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results_both{c});
    ps(c,:) = stats.pValue';
    results_both{c}
    stats.pValue
    w

    % glm with RU only
    % do model comparison
    results_RU{c} = fitglme(tbl,formula_RU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp{c} = compare(results_RU{c}, results_both{c}); % order is important -- see docs
    comp{c}
    p_comp(c,:) = comp{c}.pValue(2);
    BIC(c,:) = comp{c}.BIC';

    % glm with actRU only
    % do second model comparison
    results_actRU{c} = fitglme(tbl,formula_actRU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp2{c} = compare(results_actRU{c}, results_both{c}); % order is important -- see docs
    comp2{c}
    p_comp2(c,:) = comp2{c}.pValue(2);
    BIC2(c,:) = comp2{c}.BIC';


    % sanity check -- activations should correlate with regressor
    RU = table2array(tbl(:,'RU'));
    [r,p] = corr(RU(~exclude), act(~exclude));

    pears_rs(c,:) = r;
    pears_ps(c,:) = p;
end


save(filename, '-v7.3');


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_RU = BIC(:,1);
BIC_both = BIC(:,2);
BIC_actRU = BIC2(:,1);
table(masknames', p_uncorr, p_corr, pears_rs, pears_ps, BIC_RU, BIC_both, p_comp, BIC_actRU, p_comp2)

