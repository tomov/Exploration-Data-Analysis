% see if Badre 2012 RLPFC ROI tracks RU, option 1
% for each subject, average RU betas in ROI, then t-test across subjects
%

clear all;

data = load_data;
EXPT = exploration_expt();
glmodel = 21;
regressor = 'RU';

filename = ['badre_2012_roi_analysis_1_glm', num2str(glmodel), '.mat'];
disp(filename);

masks = badre_2012_create_masks(false);

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
    stat{i} = stats;
end

p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);

save(filename);

table(masknames', p_uncorr, p_corr, m);
