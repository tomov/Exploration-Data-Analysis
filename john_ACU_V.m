% correlate V with NAcc activity
% follow-up on john_roi_analysis
%


load('john_roi_2.mat', 'roi');
data = load_data;

load('fit_ACU_25nstarts_mixed.mat', 'results');

[tbl, lats] = data2table_ACU(roi, data, results);

formula = 'NAC ~ -1 + V';

% ignore NaN (e.g. non-existant betas for bad runs) by default
exclude = isnan(table2array(tbl(:,region)));

res = fitlme(tbl,formula,'Exclude',exclude);
[w, names, stats] = fixedEffects(res);

disp(formula);
p = stats.pValue
w
