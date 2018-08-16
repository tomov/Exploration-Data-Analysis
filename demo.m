% demo of fixedEffects and randomEffects

data = load_data;
tbl = data2table(data);


formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
%glme = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace')


[fw, fnames] = fixedEffects(glme);
[rw, rnames] = randomEffects(glme);

disp(fnames);

s = 13; % subject

% mixed effect = fixed effect + random effect for subject s
w = fw + rw((s - 1) * 3 + 1 : s * 3);

% get V,RU,VTU for subject s
tbl_s = tbl(table2array(tbl(:,'S')) == s, :);
RU = table2array(tbl_s(:,'RU'));
VTU = table2array(tbl_s(:,'VTU'));
V = table2array(tbl_s(:,'V'));

% compute decision value and prediction manually 
DV = w(1) * RU + w(2) * VTU + w(3) * V; % notice order is different -- see fnames and rnames
pred = normcdf(DV); 

% compute prediction the "right" way, make sure it's the same
y = predict(glme, tbl_s);
assert(immse(y, pred) < 1e-10);

% fit subject only
glme_s = fitglme(tbl_s,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace')

[fw_s, fnames_s] = fixedEffects(glme_s);

DV_f = fw_s(1) * RU + fw_s(2) * VTU + fw_s(3) * V;
pred_f = normcdf(DV_f);

assert(immse(y, pred_f) < 1e-10); % ! not equal !
