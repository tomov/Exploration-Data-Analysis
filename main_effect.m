function main_effect(roi_glmodel, roi_contrast, glmodel, regressor, clusterFWEcorrect, extent, Num, sphere)


printcode;

EXPT = exploration_expt();

data = load_data;


if ~exist('clusterFWEcorrect', 'var')
    clusterFWEcorrect = true;
end
if ~exist('extent', 'var')
    extent = [];
end
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end
if ~exist('sphere', 'var')
    sphere = 10; % 10 mm sphere by default
end


filename = sprintf('main_effect_roiglm%d_%s_glm%d_%s_corr=%d_extent=%d_Num=%d_s=%.1f.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, clusterFWEcorrect, extent, Num, sphere);
disp(filename);

% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num, sphere);



for i = 1:length(masks)
    mask = masks{i};
    [~, masknames{i}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, ['x', regressor], mask, s));
    end

    [h, p, ci, stats] = ttest(b);
    t = stats.tstat;
    ps(i,:) = p;
    m(i,:) = mean(b);
    disp(mask);
    p
    t
    b
    bs{i} = b;
    cis{i} = ci;
    stat{i} = stats;
end

p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);

save(filename, '-v7.3');

if exist('region', 'var')
    masknames = region';
end
table(masknames', p_uncorr, p_corr, m)
