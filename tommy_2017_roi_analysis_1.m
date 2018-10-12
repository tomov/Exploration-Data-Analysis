% see if Blanchard 2017 ROIs track TU, option 1
% for each subject, average RU betas in ROI, then t-test across subjects
%

clear all;

data = load_data;
EXPT = exploration_expt();
glmodel = 21;
regressor = 'TU';

filename = ['tommy_2017_roi_analysis_1_glm', num2str(glmodel), '.mat'];
disp(filename);

masks = tommy_2017_create_masks(false);

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

table(masknames', p_uncorr, p_corr, m)
