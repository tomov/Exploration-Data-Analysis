function cross_subject(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs, take_peak)

    % correlate neural betas with behaviora weights (betas)

% TODO dedupe with residuals_analysis.m
% TODO support what = voxel, not just sphere

if ~exist('standardize', 'var')
    standardize = false;
end
if ~exist('clusterFWEcorrect', 'var')
    clusterFWEcorrect = true;
end
if ~exist('extent', 'var')
    extent = [];
end
if ~exist('odd_runs', 'var')
    odd_runs = false;  % all runs
end
if ~exist('take_peak', 'var')
    take_peak = false; % take peak voxel? (with largest |beta|)
end

what = 'sphere';

EXPT = exploration_expt();

data = load_data;

filename = sprintf('cross_subject_roiglm%d_%s_glm%d_%s_%s_standardize=%d_corr=%d_extent=%d_tp=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, what, standardize, clusterFWEcorrect, extent, take_peak);
disp(filename);


% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent);


if odd_runs
    % behavioral weights from odd runs only
    if standardize == 1
        load results_glme_fig3_odd.mat;
    elseif standardize == 2
        load results_glme_fig3_odd_norm.mat;
    else
        load results_glme_fig3_odd_nozscore.mat;
    end
else
    % behavioral weights from all runs
    if standardize == 1
        load results_glme_fig3.mat;
    elseif standardize == 2
        load results_glme_fig3_norm.mat;
    else
        load results_glme_fig3_nozscore.mat;
    end
end

w = getEffects(results_VTURU, false);

for c = 1:length(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        if take_peak
            % take largest beta (by absolute value)
            betas = ccnl_get_beta(EXPT, glmodel, regressor, mask, s);
            [~, idx] = max(abs(betas));
            b(s) = betas(idx);
        else
            % just average betas across voxels
            b(s) = mean(ccnl_get_beta(EXPT, glmodel, regressor, mask, s));
        end
    end
    all_b{c} = b;


    if startsWith(regressor, 'RU')
        [r, p] = corr(w(:,2), b');
    elseif startsWith(regressor, 'TU')
        [r, p] = corr(w(:,3), b');
    elseif startsWith(regressor, 'V')
        [r, p] = corr(w(:,1), b');
    else
        assert(false);
    end

    disp(mask);
    r
    p

    ps(c,:) = p;
    rs(c,:) = r;
end

p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);
r = rs;

save(filename, '-v7.3');

if exist('region', 'var')
    masknames = region';
end
table(masknames', p_uncorr, p_corr, r)
