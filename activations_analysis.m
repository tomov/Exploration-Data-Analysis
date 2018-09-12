% sanity check for residuals_analysis.m -- use activations instead of residuals, and omit RU / TU -> coefficients should definitely be significant then (by definition...)
%

function [ps, results, data, region, mni, cor] = activations_analysis(glmodel, regressor, contrast, load_first_half)

outfile = ['activations_analysis_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '.mat'];

outfile

if ~exist('load_first_half', 'var')
    load_first_half = false;
end


if load_first_half
    % optionally load pre-computed activations (b/c the 2nd half doesn't work on the stupid cluster...)
    load(outfile);
else

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
            formula = 'C ~ -1 + V + VTU + actRU + (-1 + V + RU + VTU + actRU|S)'; %  omit RU
        case 'TU'
            formula = 'C ~ -1 + V + RU + VactTU + (-1 + V + RU + VTU + VactTU|S)'; %  omit VTU
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
                r = find(data(s).run(i) == runs);
                act_idx = [act_idx; idx + nTRs * (r - 1)];
            end
        end
        data(s).trial_onset_act_idx = act_idx;
    end

    % extract activations for each cluster
    %
    V_all = [];
    for s = 1:length(data) 

        act = ccnl_get_activations(EXPT, glmodel, cor, s);

        data(s).act = nan(length(data(s).run), length(region));
        [V, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
        V_all = [V_all; V(~data(s).timeout)];

        for c = 1:length(region)

            % TODO for sphere, it's
            % TODO port from context to ccnl-fmri
            % r = 10 / 1.5;
            % [~, voxels] = ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r);
            % pass voxels instead of cor(c,:)
            % then act = mean(act,2)
            %act = ccnl_get_activations(EXPT, glmodel, cor(c,:), s);

            % not all runs were used in the GLMs
            which_act = data(s).trial_onset_act_idx(~data(s).exclude); % trial onset activations
            data(s).act(~data(s).exclude,c) = squeeze(act(which_act,c,1));

            % adjust for fact that the regactsor was |RU|
            data(s).act(:,c) = data(s).act(:,c) .* (RU >= 0) + (-data(s).act(:,c)) .* (RU < 0);
        end
    end

    % in case the shit below fails
    save(outfile);
end

% fit behavioral GLM with activations
%
ps = [];
for c = 1:numel(region)
    act = [];
    exclude = logical([]);
    for s = 1:length(data)
        act = [act; data(s).act(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        exclude = [exclude; data(s).exclude(~data(s).timeout)]; % bad runs are also out (their activations are NaNs)
    end
    assert(all(isnan(act(exclude))));
    assert(all(~isnan(act(~exclude))));

    tbl = data2table(data,0,1); % exclude timeouts for fitting
    switch regressor
        case 'RU'
            actRU = act;
            tbl = [tbl table(actRU)];
        case 'TU'
            VactTU = V_all ./ act;
            tbl = [tbl table(VactTU)];
        otherwise
            assert(false);
    end

    results{c} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results{c});
    ps(c,:) = stats.pValue';
    stats.pValue
end


save(outfile);


p_uncorr = ps(:,3);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
table(region, extent, stat, mni, p_uncorr, p_corr)
