% correlate BG activity with model G and N
% uses the nonsmooth betas from GLM 23 and manually
% run GLM with G and N (fixed effects across subjects)
%

clear all;

pass = 3;

switch pass
    case 1
        % first pass -- exploratory
        masks = {'masks/striatum.nii', 'masks/pallidum.nii'}; % john_roi_1 -- first pass
        filename = 'john_roi_1.mat';
    case 2
        % second pass -- focus on ACU_mixed
        masks = {'masks/Ca.nii', 'masks/Pu.nii', 'masks/NAC.nii', 'masks/GPe.nii', 'masks/GPi.nii', 'masks/SNr.nii', 'masks/STH.nii', 'masks/SNc.nii', 'masks/VTA.nii'}; % john_roi_2 -- subcortical nuclei
        filename = 'john_roi_2.mat';
    case 3
        % third pass -- control ROIs
        masks = {'masks/v1.nii', 'masks/s1.nii', 'masks/m1.nii', 'masks/hippocampus.nii'}; % john_roi_3 -- control ROIs
        filename = 'john_roi_3.mat';
    otherwise
        assert(false);
end

roi = extract_roi_betas(masks, 'trial_onset');
save(filename, '-v7.3');

load(filename, 'roi');
data = load_data;

fitfiles = {'fit_AU_25nstarts_fixed.mat', 'fit_ACU_25nstarts_fixed.mat', 'fit_OpAL_25nstarts_fixed.mat', ...
        'fit_AU_25nstarts_mixed.mat', 'fit_ACU_25nstarts_mixed.mat', 'fit_OpAL_25nstarts_mixed.mat', ...
        'fit_AU_25nstarts_random.mat', 'fit_ACU_25nstarts_random.mat', 'fit_OpAL_25nstarts_random.mat'};
data2table_fns = {@data2table_AU, @data2table_ACU, @data2table_OpAL, ...
                @data2table_AU, @data2table_ACU, @data2table_OpAL, ...
                @data2table_AU, @data2table_ACU, @data2table_OpAL};

switch pass
    case 1
        n = length(fitfiles) * length(roi) * 2; % we're looking at all models & fits
    case 2
        n = length(roi) * 2; % we're only looking at 1 model and fit
    case 3
        n = length(roi) * 2;  % we're only looking at 1 model and fit
end

for i = 1:length(fitfiles)
    fitfile = fitfiles{i};

    load(fitfile, 'results');

    [tbl, lats] = data2table_fns{i}(roi, data, results);

    ps = []; % for G and N weights
    ws = [];
    for roi_idx = 1:length(roi)
        region = roi(roi_idx).name;
        formula = [region, ' ~ -1 + G + N'];

        % ignore NaN (e.g. non-existant betas for bad runs) by default
        exclude = isnan(table2array(tbl(:,region)));

        %roi(roi_idx).res = fitglme(tbl,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace', 'CovariancePattern','diagonal','Exclude',exclude);
        roi(roi_idx).res = fitlme(tbl,formula,'Exclude',exclude);

        [w, names, stats] = fixedEffects(roi(roi_idx).res);
        ps(roi_idx,:) = stats.pValue';
        ws(roi_idx,:) = w';
    end


    ps_corr = 1 - (1 - ps) .^ n;

    disp(fitfile);
    disp(n);
    tbl = table({roi.name}', ps(:,1), ps(:,2), ps_corr(:,1), ps_corr(:,2), ws(:,1), ws(:,2), 'VariableNames', {'ROI', 'G_uncorr', 'N_uncorr', 'G_corr', 'N_corr', 'w_G', 'w_N'});
    %tbl = table({roi.name}', ps(:,1), ps(:,2), ps(:,3), ps_corr(:,1), ps_corr(:,2), ps_corr(:,3), 'VariableNames', {'ROI', 'intercept_uncorr', 'G_uncorr', 'N_uncorr', 'intercept_corr', 'G_corr', 'N_corr'});
    disp(tbl);
end
