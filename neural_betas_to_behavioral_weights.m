function neural_betas_to_behavioral_weights(glmodel, regressor, contrast)

% TODO dedupe with residuals_analysis.m
% TODO support what = voxel, not just sphere

what = 'sphere';

outfile = ['b_to_w_', replace(contrast, ' ', '_'), '_glm', num2str(glmodel), '_', what, '.mat'];
outfile

% group-level settings
p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
direct = '+';

EXPT = exploration_expt();

data = load_data;


[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num);

radius = 10 / 1.5; % 10 mm

load results_glme_fig3_nozscore.mat;

w = getEffects(results_VTURU, false);

for c = 1:length(region)
    [~, vox] = ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), radius);

    % intersect with cluster
    c_vox = [];
    for i = 1:size(vox, 1)
        if CI(vox(i,1), vox(i,2), vox(i,3)) == CI(cor(c,1), cor(c,2), cor(c,3)) % note cluster idx != c
            c_vox = [c_vox; vox(i,:)];
        end
    end

    clear b;
    for s = 1:length(data)
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, regressor, c_vox, s));
    end


    switch regressor
        case 'RU'
            [r, p] = corr(w(:,2), b');

        case 'TU'
            [r, p] = corr(w(:,3), b');

        otherwise
            assert(false);
    end

    disp(region{c});
    r
    p

    ps(c,:) = p;
    rs(c,:) = r;
end

p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);
r = rs;

save(outfile);

table(region, p_uncorr, p_corr, r);
