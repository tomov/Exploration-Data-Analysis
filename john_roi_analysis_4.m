% See how many (unsmoothed) voxels in Str are "tuned" to G
% kind of how they do it all the time in multiunit recordings:
% find all voxels that have p < 0.01 for correlating with G
%

printcode;

clear all;

%masks = {'masks/striatum.nii', 'masks/putamen.nii', 'masks/caudate.nii', 'masks/pallidum.nii', 'masks/v1.nii', 'masks/s1.nii', 'masks/m1.nii', 'masks/hippocampus.nii'};
masks = {'masks/striatum.nii', 'masks/pallidum.nii'};

glmodels = 1019:1027; % use nosmooth betas -- nothing shows up for smooth
data = load_data;
EXPT_nosmooth = exploration_expt_nosmooth();

n = length(masks) * length(glmodels);

[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();

alpha = 0.01; % significance level

for glmodel = glmodels
    for i = 1:length(masks)
        mask = masks{i};
        [~, masknames{i}, ~] = fileparts(mask);
        m = load_mask(mask);

        clear t_G;
        clear p_G;
        for s = 1:length(data)
            t_G(s,:) = ccnl_get_tmap(EXPT_nosmooth, glmodel, 'G', mask, s);

            runs = find(goodRuns{s});
            data(s).exclude = ~ismember(data(s).run, runs) | data(s).timeout; % exclude bad runs and timeout trials
            n = sum(~data(s).exclude);
            p_G(s,:) = 2 * tcdf(abs(t_G(s,:)), n - 1, 'upper');
        end

        pos_Gs = sum(p_G < alpha, 2);
        pos_G(i,:) = mean(pos_Gs);
        pos_G_sem(i,:) = std(pos_Gs) / sqrt(size(t_G, 2));
        sz(i,:) = sum(m(:));

        fprintf('%s: %.4f +- %.4f (out of %d)\n', masknames{i}, pos_G(i), pos_G_sem(i), sz(i));
    end

    disp(['GLM ', num2str(glmodel)]);
    table(masknames', pos_G, pos_G_sem, sz, pos_G ./ sz)
end
