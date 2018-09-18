% attempt to do it just like John's ROI analysis -- from unsmoothed betas
% TODO finish it

clear all;
masks = badre_2012_create_masks(false);

%badre_roi = extract_roi_betas(masks, 'trial_onset');
%save('badre_roi.mat');

load('badre_roi.mat');
data = load_data;

tbl = data2table(data, 0, 1);

ps = []; % for G and N weights
ws = [];
for roi_idx = 1:length(roi)
	region = roi(roi_idx).name;

    b = [];
    for s = 1:length(data)
        b = [b; nanmean(roi(roi_idx).subj(s).betas, 2)];
    end

    tbl = addvars(tbl, ''); TODO fi

	formula = [region, ' ~ -1 + G + N'];

	% ignore NaN (e.g. non-existant betas for bad runs) by default
	exclude = isnan(table2array(tbl(:,region)));

	%roi(roi_idx).res = fitglme(tbl,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace', 'CovariancePattern','diagonal','Exclude',exclude);
	roi(roi_idx).res = fitlme(tbl,formula,'Exclude',exclude);

	[w, names, stats] = fixedEffects(roi(roi_idx).res);
	ps(roi_idx,:) = stats.pValue';
	ws(roi_idx,:) = w';
end

