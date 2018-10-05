% multilinear regression analysis for RU for Badre 2012 RLPFC ROI
% try to decode |RU| from multivariate ROI activity and see if it predicts
% choices better than RU from model
%
% TODO dedupe with badre_2012_activations_analysis.m


printcode;

clear all;

EXPT = exploration_expt();

data = load_data;

formula_both = 'C ~ -1 + V + RU + VTU + decRU';
formula_RU = 'C ~ -1 + V + RU + VTU';
formula_decRU = 'C ~ -1 + V + decRU + VTU';

filename = ['badre_2012_multilinear_analysis.mat'];
disp(filename);

% clusters = masks from paper
masks = badre_2012_create_masks(false);
masks = masks(1); % TODO all masks

% extract trial_onset (raw, unsmoothed) betas
%{
roi = extract_roi_betas(masks, 'trial_onset');
save(filename, '-v7.3');

load(filename, 'roi'); 

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

% clean up betas
%
for c = 1:length(roi)
    for s = 1:length(data)
        B = roi(c).subj(s).betas;
        runs = find(goodRuns{s});
        data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials
        which_nan = any(isnan(B(~data(s).exclude, :)), 1); % exclude nan voxels (ignoring bad runs and timeouts; we exclude those in the GLMs)
        B(:, which_nan) = [];
        data(s).betas{c} = B;
    end
end

% extract regressors
%
for s = 1:length(data)
    which_all = logical(ones(length(data(s).run), 1));
    [~, absRU] =  get_latents(data, s, which_all, 'abs');
    [~, RU] = get_latents(data, s, which_all, 'left');
    data(s).absRU = absRU;
    data(s).RU = RU; % for sign-correction
end

save(filename, '-v7.3');
%}

load(filename); 

for c = 1:numel(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);
    disp(mask);

    decRU = [];
    exclude = [];
    for s = 1:length(data)
        exclude = [exclude; data(s).exclude];
        X = data(s).betas{c};
        y = data(s).absRU;
        mdl = fitlm(X, y, 'exclude', data(s).exclude, 'Intercept', true);
        pred = predict(mdl, X);
        pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
        decRU = [decRU; pred];
    end
    exclude = logical(exclude);

    tbl = data2table(data, 0, 0); % include all trials; we exclude bad runs and timeouts manually
    tbl = [tbl table(decRU)];

    
    % glm with both RU and decRU
    results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results_both{c});
    ps(c,:) = stats.pValue';
    results_both{c}
    stats.pValue
    w
    names

    % glm with RU only
    % do model comparison
    results_RU{c} = fitglme(tbl,formula_RU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp{c} = compare(results_RU{c}, results_both{c}); % order is important -- see docs
    comp{c}
    p_comp(c,:) = comp{c}.pValue(2);
    BIC(c,:) = comp{c}.BIC';

    % glm with decRU only
    % do second model comparison
    results_decRU{c} = fitglme(tbl,formula_decRU,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    comp2{c} = compare(results_decRU{c}, results_both{c}); % order is important -- see docs
    comp2{c}
    p_comp2(c,:) = comp2{c}.pValue(2);
    BIC2(c,:) = comp2{c}.BIC';
end

save(filename, '-v7.3');

p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_RU = BIC(:,1);
BIC_both = BIC(:,2);
BIC_decRU = BIC2(:,1);
table(masknames', p_uncorr, p_corr, BIC_RU, BIC_both, p_comp, BIC_decRU, p_comp2)

