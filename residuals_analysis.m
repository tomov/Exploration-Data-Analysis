function [ps, results, data, region, mni, cor] = residuals_analysis(glmodel, regressor, contrast)


% group-level settings
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
direct = '+';

EXPT = exploration_expt();

data = load_data;

switch regressor
    case 'RU'
        formula = 'C ~ -1 + V + RU + VTU + resRU + (-1 + V + RU + VTU + resRU|S)';
    case 'TU'
        formula = 'C ~ -1 + V + RU + VTU + VresTU + (-1 + V + RU + VTU + VresTU|S)';
    otherwise
        assert(false);
end

[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);

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

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

% extract residuals for each cluster
%
V_all = [];
for s = 1:length(data) 

    res = ccnl_get_residuals(EXPT, glmodel, cor, s);

    data(s).res = nan(length(data(s).run), length(region));
    [V, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
    V_all = [V_all; V(~data(s).timeout)];

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
        data(s).res(which_trials,c) = squeeze(res(which_res,c,1));

        % adjust for fact that the regressor was |RU|
        data(s).res(:,c) = data(s).res(:,c) .* (RU >= 0) + (-data(s).res(:,c)) .* (RU < 0);
    end
end

% in case the shit below fails
save(['residuals_analysis_', regressor, '_glm', num2str(glmodel), '.mat']);

% fit behavioral GLM with residuals
%
ps = [];
for c = 1:numel(region)
    res = [];
    for s = 1:length(data)
        res = [res; data(s).res(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
    end

    tbl = data2table(data,0,1); % exclude timeouts for fitting
    switch regressor
        case 'RU'
            resRU = res;
            tbl = [tbl table(resRU)];
        case 'TU'
            VresTU = V_all ./ res;
            tbl = [tbl table(VresTU)];
        otherwise
            assert(false);
    end


    results{c} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal');
    [w, names, stats] = fixedEffects(results{c});
    ps(c) = stats(4).pValue;
end


save(['residuals_analysis_', regressor, '_glm', num2str(glmodel), '.mat']);
