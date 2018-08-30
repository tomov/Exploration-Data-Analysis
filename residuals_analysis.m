%function [ps, results, data] = residuals_analysis()


% group-level settings
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
direct = '+';

EXPT = exploration_expt();
glmodel = 10;

data = load_data;

% TODO for TU, it's TU - trial
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, 'RU', p, direct, alpha, Dis, Num);

% find peak of HRF
hrf = spm_hrf(0.001);
[~,hrf_offset] = max(hrf);
hrf_offset = hrf_offset / 1000;

nTRs = 242;
TR = EXPT.TR;
trs = TR/2 : TR : nTRs * TR;

% find closest TR to each trial onset (adjusted for HRF f'n)
for s = 1:length(data)
    res_idx = [];
    for i = 1:length(data(s).trial_onset)
        [~, idx] = min(abs(trs - (data(s).trial_onset(i) + hrf_offset)));
        res_idx = [res_idx; idx];
    end
    data(s).trial_onset_res_idx = res_idx;
end

% TODO for TU it's VresTU
formula = 'C ~ -1 + V + RU + VTU + resRU + (-1 + V + RU + VTU + resRU|S)';

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

% extract residuals for each cluster
%
ps = [];
for s = 1:length(data) 

    res = ccnl_get_residuals(EXPT, glmodel, cor, s);

    data(s).resRU = nan(length(data(s).run), length(region));
    [~, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');

    for c = 1:length(region)

        % TODO for sphere, it's
        % TODO port from context to ccnl-fmri
        % r = 10 / 1.5;
        % [~, voxels] = ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r);
        % pass voxels instead of cor(c,:)
        % then res = mean(res,2)
        %res = ccnl_get_residuals(EXPT, glmodel, cor(c,:), s);

        % not all runs were used in the GLMs
        which_trials = ismember(data(s).run, find(goodRuns{s}));
        which_res = data(s).trial_onset_res_idx(which_trials);
        data(s).resRU(which_trials,c) = squeeze(res(which_res,c,1));

        % adjust for fact that the regressor was |RU|
        data(s).resRU(:,c) = data(s).resRU(:,c) .* (RU >= 0) + (-data(s).resRU(:,c)) .* (RU < 0);
    end
end

% fit behavioral GLM with residuals
%
for c = 1:numel(region)
    resRU = [];
    for s = 1:length(data)
        resRU = [resRU; data(s).resRU(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
    end

    tbl = data2table(data,0,1); % exclude timeouts for fitting
    % TODO for TU, it's VresTU = table2array(tbl('V')) ./ resTU
    tbl = [tbl table(resRU)];

    results(c) = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal');
    [w, names, stats] = fixedEffects(results);
    assert(strcmp(names{4}, 'resRU'));
    ps(c) = stats(4).pValue;
end

save('results_analysis.mat', 'ps', 'results', 'data');
