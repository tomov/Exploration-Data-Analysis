% see if Badre 2012 RLPFC ROI tracks RU, option 2
% for each voxel in ROI, t-test RU beta across subjects
% then Bonferonni correct
%

clear all;

data = load_data;
EXPT = exploration_expt();
glmodel = 21;
regressor = 'RU';

alpha = 0.05;

masks = badre_2012_create_masks(false);

for i = 1:length(masks)
    mask = masks{i};
    [~, masknames{i}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        b(s, :) = ccnl_get_beta(EXPT, glmodel, regressor, mask, s);
    end

    [h, p, ci, stats] = ttest(b);

    p = 1 - (1 - p) .^ length(p);
    t = stats.tstat;
    disp(mask);
    disp(any(p < alpha));
    p
    t
    mean(b, 1)

    all_p{i} = p;
    all_stat{i} = stats;
    min_p(i,:) = min(p);
end


save('badre_2012_roi_analysis_2.mat');

table(masknames', min_p);
