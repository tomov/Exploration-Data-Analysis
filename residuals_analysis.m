% see if residuals for peak RU / TU voxels contribute to choice prediction beyond model estimates
%

function [ps, results, data, region, mni, cor] = residuals_analysis(glmodel, regressor, contrast, what, load_first_half)

printcode;

if ~exist('load_first_half', 'var')
    load_first_half = false;
end
if ~exist('what', 'var')
    what = 'voxel';
end

outfile = ['residuals_analysis_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '_', what, '.mat'];
outfile



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
            formula = 'C ~ -1 + V + RU + VTU + resRU';
        case 'TU'
            formula = 'C ~ -1 + V + RU + VTU + VresTU';
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

        switch what
            case 'voxel'
                res = ccnl_get_residuals(EXPT, glmodel, cor, s);

            case 'sphere'
                r = 10 / 1.5; % 10 mm radius
                clear res;
                for c = 1:length(region)
                    [~, vox] = ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r);

                    % intersect with cluster
                    c_vox = [];
                    for i = 1:size(vox, 1)
                        if CI(vox(i,1), vox(i,2), vox(i,3)) == CI(cor(c,1), cor(c,2), cor(c,3)) % note cluster idx != c
                            c_vox = [c_vox; vox(i,:)];
                        end
                    end
                    res(:,c) = mean(ccnl_get_residuals(EXPT, glmodel, c_vox, s), 2);
                end

            otherwise
                assert(false, 'what should be voxel or sphere');
        end

        data(s).res = nan(length(data(s).run), length(region));
        [V, RU] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');
        V_all = [V_all; V(~data(s).timeout)];

        for c = 1:length(region)
            % not all runs were used in the GLMs
            which_res = data(s).trial_onset_res_idx(~data(s).exclude); % trial onset residuals
            data(s).res(~data(s).exclude,c) = squeeze(res(which_res,c,1));

            % adjust for fact that the regressor was |RU| 
            if strcmp(regressor, 'RU')
                data(s).res(:,c) = data(s).res(:,c) .* (RU >= 0) + (-data(s).res(:,c)) .* (RU < 0);
            end
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

    % sanity check -- residuals should NOT correlate with regressor
    switch regressor
        case 'RU'
            RU = table2array(tbl(:,'RU'));
            [r,p] = corr(RU(~exclude), res(~exclude));
        case 'TU'
            TU = table2array(tbl(:,'TU'));
            [r,p] = corr(TU(~exclude), res(~exclude));
        otherwise
            assert(false);
    end
    pears_rs(c,:) = r;
    pears_ps(c,:) = p;
end


save(outfile);


p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
table(region, extent, stat, mni, p_uncorr, p_corr, pears_rs, pears_ps)
