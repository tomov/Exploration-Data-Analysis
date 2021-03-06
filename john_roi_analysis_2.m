% correlate BG with model G and N
% average (non)smooth betas from GLMs 1001..1009 within the ROI
% do a t-test across subjects
% similar to badre_2012_roi_analysis_1.m
%

printcode;

clear all;

%masks = {'masks/striatum.nii', 'masks/putamen.nii', 'masks/caudate.nii', 'masks/pallidum.nii', 'masks/v1.nii', 'masks/s1.nii', 'masks/m1.nii', 'masks/hippocampus.nii'};
masks = {'masks/striatum.nii', 'masks/pallidum.nii'};

glmodels = 1001:1009; % use nosmooth betas -- nothing shows up for smooth
data = load_data;
EXPT = exploration_expt();

n = length(masks) * length(glmodels);

for glmodel = glmodels
    for i = 1:length(masks)
        mask = masks{i};
        [~, masknames{i}, ~] = fileparts(mask);

        clear b_G;
        clear b_N;
        for s = 1:length(data)
            b_G(s) = nanmean(ccnl_get_beta_nosmooth(EXPT, glmodel, 'G', mask, s));
            b_N(s) = nanmean(ccnl_get_beta_nosmooth(EXPT, glmodel, 'N', mask, s));
        end

        [h, p, ci, stats] = ttest(b_G);
        t_G(i,:) = stats.tstat;
        ps_G(i,:) = p;
        w_G(i,:) = mean(b_G);

        [h, p, ci, stats] = ttest(b_N);
        t_N(i,:) = stats.tstat;
        ps_N(i,:) = p;
        w_N(i,:) = mean(b_N);

        disp([mask, ', G']);
        t_G(i,:)
        p
        w_G(i,:)

        disp([mask, ', N']);
        t_N(i,:)
        p
        w_N(i,:)
    end

    G_uncorr = ps_G;
    G_corr = 1 - (1 - ps_G) .^ length(ps_G);

    N_uncorr = ps_N;
    N_corr = 1 - (1 - ps_N) .^ length(ps_N);

    disp(['GLM ', num2str(glmodel)]);
    table(masknames', G_uncorr, N_uncorr, G_corr, N_corr, w_G, w_N, t_G, t_N)
end
