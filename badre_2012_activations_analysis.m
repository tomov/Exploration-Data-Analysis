% activations analysis for RU for Badre 2012 RLPFC ROI
% TODO dedupe with activations_analysis.m
% TODO dedupe with badre_2012_residuals_analysis_glm.m

function badre_2012_activations_analysis(glmodel, normalize)

printcode;

EXPT = exploration_expt();

if ~exist('glmodel', 'var')
    glmodel = 21;
end
if ~exist('normalize', 'var')
    normalize = true; % divide each activation by the corresponding beta
end

data = load_data;

formula = 'C ~ -1 + V + RU + VTU + actRU';
formula_RU = 'C ~ -1 + V + RU + VTU';

if normalize
    filename = ['badre_2012_activations_analysis_glm', num2str(glmodel), '_normalized.mat'];
else
    filename = ['badre_2012_activations_analysis_glm', num2str(glmodel), '.mat'];
end
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

% get betas to (optionally) normalize activations in each run
for c = 1:length(masks)
    mask = masks{c};
    m = load_mask(mask);
    cnt = sum(m(:));

    for s = 1:length(data)
        runs = find(goodRuns{s});
        data(s).b{c} = nan(length(data(s).run), cnt);

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
    end
end

% extract activations for each cluster
%
V_all = [];
for s = 1:length(data)

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
        % not all runs were used in the GLMs
        which_act = data(s).trial_onset_act_idx(~data(s).exclude); % trial onset activations
        act{c} = act{c}(which_act,:); % only consider 1 activation for each trial

        if normalize % alternatively, do act - b0, where b0 is averaged across runs (i.e. for entire subject)
            act{c} = (act{c} - data(s).b0{c}(~data(s).exclude)) ./ data(s).b{c}(~data(s).exclude);
        end

        data(s).act(~data(s).exclude,c) = mean(act{c}, 2);

        % adjust for fact that the regactsor was |RU|
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

    results{c} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results{c});
    ps(c,:) = stats.pValue';
    results{c}
    stats.pValue
    w

    % model comparison with original formula
    results_RU{c} = fitglme(tbl,formula_RU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp{c} = compare(results_RU{c}, results{c}); % order is important -- see docs
    comp{c}
    p_comp(c,:) = comp{c}.pValue(2);
    BIC(c,:) = comp{c}.BIC';


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
BIC_actRU = BIC(:,2);
table(masknames', p_uncorr, p_corr, pears_rs, pears_ps, BIC_RU, BIC_actRU, p_comp)

