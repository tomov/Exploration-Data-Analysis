% correlate BG activity with model G and N
%

%load('john_roi.mat', 'roi');
data = load_data;
load('fit_AU_25nstarts_fixed.mat', 'results');

tbl = data2table_AU(roi, data, results);

for roi_idx = 1:length(roi)
    formula = [roi(roi_idx).name, ' ~ 1 + G + N + (1 + G + N | S)'];

    % fitglme ignores NaN (e.g. non-existant betas for bad runs) by default
    roi(roi_idx).res = fitglme(tbl,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace', 'CovariancePattern','diagonal');
end
