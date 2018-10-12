% correlate BG with model G and N
% use nonsmooth betas averaged within the ROI
% do a t-test across subjects
% similar to badre_2012_roi_analysis_1.m
%

printcode;

clear all;

%masks = {'masks/striatum.nii', 'masks/putamen.nii', 'masks/caudate.nii', 'masks/pallidum.nii', 'masks/v1.nii', 'masks/s1.nii', 'masks/m1.nii', 'masks/hippocampus.nii'};
masks = {'masks/striatum.nii', 'masks/pallidum.nii'};

glmodels = 1019:1027;
data = load_data;
EXPT = exploration_expt_nosmooth();

n = length(masks) * length(glmodels);

for glmodel = glmodels
    for i = 1:length(masks)
        mask = masks{i};
        [~, masknames{i}, ~] = fileparts(mask);

        b_G = ccnl_roi_glm(EXPT, glmodel, mask, 'G');

        [h, p, ci, stats] = ttest(b_G);
        t_G(i,:) = stats.tstat;
        ps_G(i,:) = p;
        w_G(i,:) = mean(b_G);

        %[h, p, ci, stats] = ttest(b_N);
        %t_N(i,:) = stats.tstat;
        %ps_N(i,:) = p;
        %w_N(i,:) = mean(b_N);

        disp([mask, ', G']);
        t_G(i,:)
        p
        w_G(i,:)

        %disp([mask, ', N']);
        %t_N(i,:)
        %p
        %w_N(i,:)
    end

    G_uncorr = ps_G;
    G_corr = 1 - (1 - ps_G) .^ length(ps_G);

    %N_uncorr = ps_N;
    %N_corr = 1 - (1 - ps_N) .^ length(ps_N);

    disp(['GLM ', num2str(glmodel)]);
    table(masknames', G_uncorr, G_corr, w_G, t_G)
    %table(masknames', G_uncorr, N_uncorr, G_corr, N_corr, w_G, w_N, t_G, t_N)
end
