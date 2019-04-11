% correlate residuals between DV and RU and TU ROIs
% copied from roi_corr.m
%

clear all;

DV_glm = 47;
VTURU_glm = 36;


data = load_data;

EXPT = exploration_expt();
clusterFWEcorrect = false;
extent = 100;
Num = 1;

% get DV ROI from GLM 47
[DV_masks, ~] = get_masks(DV_glm, 'DV', clusterFWEcorrect, extent, Num);
DV_mask = DV_masks{1};

for s = 1:length(data)
    DV_res{s} = mean(ccnl_get_residuals(EXPT, DV_glm, DV_mask, s), 2);
end


% get RU ROIs from GLM 36
[RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);

% get TU ROIs from GLM 36
[TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);



for c = 1:length(RU_masks)

    disp(RU_masks{c});

    rs = [];
    for s = 1:length(data)
        RU_res{s} = mean(ccnl_get_residuals(EXPT, VTURU_glm, RU_masks{c}, s), 2);

        [r, p] = corr(RU_res{s}, DV_res{s});
        rs(s) = r;
    end

    rs

    rs = atanh(rs);
    [h, p, ci, stats] = ttest(rs)
end
