% see if |DV| predicts reaction times



data = load_data;
tbl = data2table(data,0,1); % don't standardize! TODO never standardize anywhere

DV_all = [];
V_all = [];
timeout = [];
for s = 1:length(data)

    % get behavioral regressors
    [V, RU, TU, VTU, DV] = get_latents(data, s, logical(ones(length(data(s).run), 1)), 'left');

    V_all = [V_all; V];
    DV_all = [DV_all; DV];
    timeout = [timeout; data(s).timeout];
end
V = V_all(~timeout);
DV = DV_all(~timeout);

absDV = abs(DV);
absV = abs(V);

assert(immse(V, tbl.V) < 1e-9);

tbl = [tbl table(DV, absDV, absV)];



formula_DV = 'rt ~ 1 + absDV + (1 + absDV | S)';
result_DV = fitglme(tbl, formula_DV, 'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace');
result_DV

[beta, names, stats] = fixedEffects(result_DV);
H = [0 1];
[p, F, DF1, DF2] = coefTest(result_DV, H);
fprintf('fitglme "%s": DV beta = %f, p = %f, F(%d,%d) = %f\n', formula_DV, H * beta, p, DF1, DF2, F);





formula_V = 'rt ~ 1 + absV + (1 + absV | S)';
result_V = fitglme(tbl, formula_V, 'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace');
result_V

[beta, names, stats] = fixedEffects(result_V);
H = [0 1];
[p, F, DF1, DF2] = coefTest(result_V, H);
fprintf('fitglme "%s": V beta = %f, p = %f, F(%d,%d) = %f\n', formula_V, H * beta, p, DF1, DF2, F);







formula_VDV = 'rt ~ 1 + absDV + absV + (1 + absDV + absV | S)';
result_VDV = fitglme(tbl, formula_VDV, 'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace');
result_VDV

[beta, names, stats] = fixedEffects(result_VDV);
H = [0 1 0];
[p, F, DF1, DF2] = coefTest(result_VDV, H);
fprintf('fitglme "%s": DV beta = %f, p = %f, F(%d,%d) = %f\n', formula_VDV, H * beta, p, DF1, DF2, F);


comp = compare(result_V, result_VDV);
comp


