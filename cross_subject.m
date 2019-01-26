function cross_subject(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)

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

what = 'sphere';

EXPT = exploration_expt();

data = load_data;

filename = sprintf('cross_subject_roiglm%d_%s_glm%d_%s_%s_standardize=%d_corr=%d_extent=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, what, standardize, clusterFWEcorrect, extent);
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
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, regressor, mask, s));
    end
    all_b{c} = b;


    switch regressor
        case 'RU'
            [r, p] = corr(w(:,2), b');

        case 'TU'
            [r, p] = corr(w(:,3), b');

        otherwise
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
