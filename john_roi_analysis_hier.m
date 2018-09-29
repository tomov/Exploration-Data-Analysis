% correlate BG activity with model G and N
% same as john_roi_analysis except hierarchical fit, then correlating weights for G and N with a and b for subjects (got nothing...)
%


load('john_roi_1.mat', 'roi');
data = load_data;

fitfiles = {'fit_AU_25nstarts_fixed.mat', 'fit_ACU_25nstarts_fixed.mat', 'fit_OpAL_25nstarts_fixed.mat', ...
        'fit_AU_25nstarts_mixed.mat', 'fit_ACU_25nstarts_mixed.mat', 'fit_OpAL_25nstarts_mixed.mat', ...
        'fit_AU_25nstarts_random.mat', 'fit_ACU_25nstarts_random.mat', 'fit_OpAL_25nstarts_random.mat'};
data2table_fns = {@data2table_AU, @data2table_ACU, @data2table_OpAL, ...
                @data2table_AU, @data2table_ACU, @data2table_OpAL, ...
                @data2table_AU, @data2table_ACU, @data2table_OpAL};

n = length(fitfiles) * length(roi);

for i = 1:length(fitfiles)
    fitfile = fitfiles{i};

    load(fitfile, 'results');

    [tbl, lats] = data2table_fns{i}(roi, data, results);

    ps = []; % for G and N weights
    ras = []; % for correlation btwn G weights and a
    rbs = []; % for correlation btwn N weights and b
    pas = []; % for correlation btwn G weights and a
    pbs = []; % for correlation btwn N weights and b
    for roi_idx = 1:length(roi)
        region = roi(roi_idx).name;
        formula = [region, ' ~ -1 + G + N + (-1 + G + N | S)'];

        exclude = isnan(table2array(tbl(:,region)));

        % fitglme ignores NaN (e.g. non-existant betas for bad runs) by default
        %roi(roi_idx).res = fitglme(tbl,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace', 'CovariancePattern','diagonal','Exclude',exclude);
        roi(roi_idx).res = fitlme(tbl,formula,'Exclude',exclude);

        [w, names, stats] = fixedEffects(roi(roi_idx).res);
        ps(roi_idx,:) = stats.pValue';

        if ~isempty(strfind(fitfile, 'mixed')) || ~isempty(strfind(fitfile, 'random'))
            [w_rand, names_rand] = randomEffects(roi(roi_idx).res);
            w_G = w_rand(1:2:end) + w(1);
            w_N = w_rand(2:2:end) + w(2);
            a = [lats.a]';
            b = [lats.b]';
            [ra, pa] = corr(w_G, a);
            [rb, pb] = corr(w_N, b);
            ras = [ras; ra];
            rbs = [rbs; rb];
            pas = [pas; pa];
            pbs = [pbs; pb];
        end
    end

    ps_corr = 1 - (1 - ps) .^ n;

    disp(fitfile);
    if length(ras) > 0
        tbl = table({roi.name}', ps(:,1), ps(:,2), ps_corr(:,1), ps_corr(:,2), ras, pas, rbs, pbs, 'VariableNames', {'ROI', 'G_uncorr', 'N_uncorr', 'G_corr', 'N_corr', 'r_G', 'p_G', 'r_N', 'p_N'});
    else
        tbl = table({roi.name}', ps(:,1), ps(:,2), ps_corr(:,1), ps_corr(:,2), 'VariableNames', {'ROI', 'G_uncorr', 'N_uncorr', 'G_corr', 'N_corr'});
    end
    %tbl = table({roi.name}', ps(:,1), ps(:,2), ps(:,3), ps_corr(:,1), ps_corr(:,2), ps_corr(:,3), 'VariableNames', {'ROI', 'intercept_uncorr', 'G_uncorr', 'N_uncorr', 'intercept_corr', 'G_corr', 'N_corr'});
    disp(tbl);
end
