function cross_subject(glmodel, regressor, contrast, standardize)

% TODO dedupe with residuals_analysis.m
% TODO support what = voxel, not just sphere

if ~exist('standardize', 'var')
    standardize = false;
end

what = 'sphere';

EXPT = exploration_expt();

data = load_data;

filename = sprintf('cross_subject_glm%d_%s_%s_%s.mat', glmodel, regressor, replace(contrast, ' ', '_'), what);
disp(filename);

% get ROI masks
switch contrast
    case 'badre'
        % clusters = masks from paper
        masks = badre_2012_create_masks(false);
        %masks = masks(1); % TODO use all masks

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
end



if standardize
    load results_glme_fig3_nozscore.mat;
else
    load results_glme_fig3.mat;
end
w = getEffects(results_VTURU, false);

for c = 1:length(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);

    %{
    [~, vox] = ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), radius);

    % intersect with cluster
    c_vox = [];
    for i = 1:size(vox, 1)
        if CI(vox(i,1), vox(i,2), vox(i,3)) == CI(cor(c,1), cor(c,2), cor(c,3)) % note cluster idx != c
            c_vox = [c_vox; vox(i,:)];
        end
    end
    %}

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
