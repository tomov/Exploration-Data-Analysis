function cross_subj_vbm(glmodel, regressor, contrast, standardize, clusterFWEcorrect, extent, Num)
% correlate grey matter thickness with weight for given regressor (e.g. RU)
% -> maybe subjects who have thicker RLPFC rely on RU more
%

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
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end

what = 'sphere';

EXPT = exploration_expt();

data = load_data;

filename = sprintf('cross_subj_vbm_glm%d_%s_%s_%s_standardize=%d_corr=%d_extent=%d_Num=%d.mat', glmodel, regressor, replace(contrast, ' ', '_'), what, standardize, clusterFWEcorrect, extent, Num);
disp(filename);

% get ROIs
[masks, region] = get_masks(glmodel, contrast, clusterFWEcorrect, extent, Num);



if standardize == 1
    load results_glme_fig3.mat;
elseif standardize == 2
    load results_glme_fig3_norm.mat;
else
    load results_glme_fig3_nozscore.mat;
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

    b = ccnl_vbm(EXPT, mask);
    all_b{c} = b;


    switch regressor
        case 'RU'
            [r, p] = corr(w(:,2), b');

        case 'TU'
            [r, p] = corr(w(:,3), b');

        case 'V'
            [r, p] = corr(w(:,1), b');

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
