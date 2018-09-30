% residuals analysis for RU for Badre 2012 RLPFC ROI
% TODO dedupe with residuals_analysis.m

clear all;


EXPT = exploration_expt();

glmodel = 21;

data = load_data;

normalize = true; % divide each residual by the corresponding beta

formula = 'C ~ -1 + V + RU + VTU + resRU';

filename = ['badre_2012_residuals_analysis_glm', num2str(glmodel), '.mat'];
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
    res_idx = [];
    runs = find(goodRuns{s});
    data(s).exclude = ~ismember(data(s).run, runs); % exclude bad runs
    for i = 1:length(data(s).trial_onset)
        [~, idx] = min(abs(trs - (data(s).trial_onset(i) + hrf_offset)));
        if data(s).exclude(i)
            res_idx = [res_idx; NaN];
        else
            r = find(data(s).run(i) == runs); % scan session idx in GLM 
            res_idx = [res_idx; idx + nTRs * (r - 1)];
        end
    end
    data(s).trial_onset_res_idx = res_idx;
end


% clusters = masks from paper
masks = badre_2012_create_masks(false);

% get betas to (optionally) normalize residuals in each run
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
                reg = ['Sn(', num2str(r), ') trial_onsetxRU'];
                fprintf('  c = %s, s = %d, run = %d, r = %d, reg = %s\n', mask, s, run, r, reg);
                b = ccnl_get_beta(EXPT, glmodel, reg, mask, s);
                data(s).b{c}(data(s).run == run, :) = repmat(b, sum(data(s).run == run), 1);
            end
        end
    end
end

% extract residuals for each cluster
%
V_all = [];
for s = 1:length(data)

    clear res;
    for c = 1:length(masks)
        mask = masks{c};
        [~, masknames{c}, ~] = fileparts(mask);

        res{c} = ccnl_get_residuals(EXPT, glmodel, mask, s);
        data(s).all_res{c} = res{c};

    end

    data(s).res = nan(length(data(s).run), length(masks));
    [V, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    V_all = [V_all; V(~data(s).timeout)];

    for c = 1:length(masks)
        % not all runs were used in the GLMs
        which_res = data(s).trial_onset_res_idx(~data(s).exclude); % trial onset residuals
        res{c} = res{c}(which_res,:); % only consider 1 residual for each trial

        if normalize
            res{c} = res{c} ./ data(s).b{c}(~data(s).exclude);
        end

        data(s).res(~data(s).exclude,c) = mean(res{c}, 2);

        % adjust for fact that the regressor was |RU|
        if glmodel == 21
            data(s).res(:,c) = data(s).res(:,c) .* (RU >= 0) + (-data(s).res(:,c)) .* (RU < 0);
        end
    end
end

save(filename);


% fit behavioral GLM with residuals
%
ps = [];
for c = 1:numel(masks)
    res = [];
    exclude = logical([]);
    for s = 1:length(data)
        res = [res; data(s).res(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        exclude = [exclude; data(s).exclude(~data(s).timeout)]; % bad runs are also out (their residuals are NaNs)
    end
    assert(all(isnan(res(exclude))));
    assert(all(~isnan(res(~exclude))));

    tbl = data2table(data,0,1); % exclude timeouts for fitting
    resRU = res;
    tbl = [tbl table(resRU)];

    results{c} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results{c});
    ps(c,:) = stats.pValue';
    stats.pValue

    % sanity check -- residuals should NOT correlate with regressor
    RU = table2array(tbl(:,'RU'));
    [r,p] = corr(RU(~exclude), res(~exclude));

    pears_rs(c,:) = r;
    pears_ps(c,:) = p;
end


save(filename);


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
table(masknames', p_uncorr, p_corr, pears_rs, pears_ps)

