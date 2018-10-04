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

% extract trial_onset (raw, unsmoothed) betas
roi = extract_roi_betas(masks, 'trial_onset');
save(filename, '-v7.3');

load(filename, 'roi'); 

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

for roi_idx = 1:length(roi)
    for s = 1:length(s)
        B = roi(roi_idx).subj(s).betas;
        runs = find(goodRuns{s});
        data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials
        which_nan = any(isnan(B(~data(s).exclude, :), 1); % exclude nan voxels (ignoring bad runs and timeouts)
        B(:, which_nan) = [];
        data(s).betas{roi_idx} = B;
    end
end

for s = 1:length(s)
    which_all = logical(ones(length(data(s).run), 1));
    [~, absRU] =  get_latents(data, s, which_all, 'abs');
    [~, RU] = get_latents(data, s, which_all, 'left');
    data(s).absRU = absRU;
    data(s).RU = RU;
end

for roi_idx = 1:length(roi)
    decRU = [];
    exclude = [];
    for s = 1:length(data)
        exclude = [exclude; data(s).exclude];
        mdl = fitlm(data(s).betas{roi_idx}, data(s).absRU, 'exclude', data(s).exclude, 'Intercept', true);
        pred = predict(mdl, data(s).betas{roi_idx}); % TODO do leave-one-run-out to avoid overfitting and just recovering |RU| 
        pred = pred .* (data(s).RU >= 0) + (-pred) .* (data(s).RU < 0); % adjust for fact that we decode |RU|
        decRU = [decRU; pred];
    end

    tbl = data2table(data, 0, 0); % include all trials; we exclude bad runs and timeouts manually
    tbl = tbl(~exclude, :);
    tbl = [tbl table(decRU)];

    
    % glm with both RU and decRU
    results_both{c} = fitglme(tbl,formula_both,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal', 'Exclude',exclude);
    [w, names, stats] = fixedEffects(results_both{c});
    ps(c,:) = stats.pValue';
    results_both{c}
    stats.pValue
    w

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


    % sanity check -- activations should correlate with regressor
    RU = table2array(tbl(:,'RU'));
    [r,p] = corr(RU(~exclude), act(~exclude));

    pears_rs(c,:) = r;
    pears_ps(c,:) = p;
end

save(filename, '-v7.3');

p_uncorr = ps(:,4);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_RU = BIC(:,1);
BIC_both = BIC(:,2);
BIC_decRU = BIC2(:,1);
table(masknames', p_uncorr, p_corr, pears_rs, pears_ps, BIC_RU, BIC_both, p_comp, BIC_decRU, p_comp2)

