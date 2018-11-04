function main_effect(glmodel, regressor, contrast)

printcode;

EXPT = exploration_expt();

data = load_data;

filename = sprintf('main_effect_glm%d_%s_%s.mat', glmodel, regressor, replace(contrast, ' ', '_'));
disp(filename);


% get ROI masks
switch contrast
    case 'badre'
        % clusters = masks from paper
        masks = badre_2012_create_masks(false);
        %masks = masks(1); % TODO use all masks

    case 'dlpfc'
        % clusters = masks from paper
        masks = dlpfc_2012_create_masks(false);

    case 'tommy'
        % clusters = masks from paper
        masks = tommy_2017_create_masks(false);

    otherwise
        % group-level settings
        p = 0.001;
        alpha = 0.05;
        Dis = 20;
        Num = 1; % # peak voxels per cluster; default in bspmview is 3
        direct = '+';

        [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);

        r = 10 / 1.5; % 10 mm radius

        % create spherical masks around peak voxel of each cluster (intersected with cluster)
        %
        for c = 1:length(region)
            masks{c} = sprintf('sphere_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));
            cmask = CI == CI(cor(c,1), cor(c,2), cor(c,3));
            ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c}, cmask);
        end
        masknames = region';
end


clear stat;


for i = 1:length(masks)
    mask = masks{i};
    [~, masknames{i}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, regressor, mask, s));
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
