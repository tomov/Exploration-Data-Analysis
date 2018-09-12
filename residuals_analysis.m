% see if residuals for peak RU / TU voxels contribute to choice prediction beyond model estimates
%

function [ps, results, data, region, mni, cor] = residuals_analysis(glmodel, regressor, contrast, load_first_half)

outfile = ['residuals_analysis_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '.mat'];

outfile

if ~exist('load_first_half', 'var')
    load_first_half = false;
end


if load_first_half
    % optionally load pre-computed residuals (b/c the 2nd half doesn't work on the stupid cluster...)
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
                r = find(data(s).run(i) == runs);
                res_idx = [res_idx; idx + nTRs * (r - 1)];
            end
        end
        data(s).trial_onset_res_idx = res_idx;
    end

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
            which_res = data(s).trial_onset_res_idx(~data(s).exclude); % trial onset residuals
            data(s).res(~data(s).exclude,c) = squeeze(res(which_res,c,1));

            % adjust for fact that the regressor was |RU|
            data(s).res(:,c) = data(s).res(:,c) .* (RU >= 0) + (-data(s).res(:,c)) .* (RU < 0);
        end
    end

    % in case the shit below fails
    save(outfile);
end

% fit behavioral GLM with residuals
%
ps = [];
for c = 1:numel(region)
    res = [];
    exclude = logical([]);
    for s = 1:length(data)
        res = [res; data(s).res(~data(s).timeout, c)]; % even though neural GLMs includes timeouts, we exclude them for fitting the behavioral GLMs
        exclude = [exclude; data(s).exclude(~data(s).timeout)]; % bad runs are also out (their residuals are NaNs)
    end
    assert(all(isnan(res(exclude))));
    assert(all(~isnan(res(~exclude))));

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

    results{c} = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results{c});
    ps(c,:) = stats.pValue';
    stats.pValue
end


save(outfile);


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
table(region, extent, stat, mni, p_uncorr, p_corr)
