% correlate BG activity with model G and N
%


%load('john_roi.mat', 'roi');
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
    ws = [];
    for roi_idx = 1:length(roi)
        region = roi(roi_idx).name;
        formula = [region, ' ~ -1 + G + N'];

        exclude = isnan(table2array(tbl(:,region)));

        % fitglme ignores NaN (e.g. non-existant betas for bad runs) by default
        %roi(roi_idx).res = fitglme(tbl,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace', 'CovariancePattern','diagonal','Exclude',exclude);
        roi(roi_idx).res = fitlme(tbl,formula,'Exclude',exclude);

        [w, names, stats] = fixedEffects(roi(roi_idx).res);
        ps(roi_idx,:) = stats.pValue';
        ws(roi_idx,:) = w';
    end

    ps_corr = 1 - (1 - ps) .^ n;

    disp(fitfile);
    tbl = table({roi.name}', ps(:,1), ps(:,2), ps_corr(:,1), ps_corr(:,2), ws(:,1), ws(:,2), 'VariableNames', {'ROI', 'G_uncorr', 'N_uncorr', 'G_corr', 'N_corr', 'w_G', 'w_N'});
    %tbl = table({roi.name}', ps(:,1), ps(:,2), ps(:,3), ps_corr(:,1), ps_corr(:,2), ps_corr(:,3), 'VariableNames', {'ROI', 'intercept_uncorr', 'G_uncorr', 'N_uncorr', 'intercept_corr', 'G_corr', 'N_corr'});
    disp(tbl);
end
