% WRONG -- this is super circular;
% the RU and TU ROIs correlate with RU and TU
% the DV ROI correlates with DV which correlates with RU and TU
% so ofc their beta series will be correlated doh

% correlate beta series between DV and RU and TU ROIs
%

clear all;

DV_glm = 47;
VTURU_glm = 36;
beta_series_glm = 23;


data = load_data;

EXPT = exploration_expt();
clusterFWEcorrect = false;
extent = 100;
Num = 1;

% get DV ROI from GLM 47
[DV_masks, ~] = get_masks(DV_glm, 'DV', clusterFWEcorrect, extent, Num);
DV_mask = DV_masks{1};

for s = 1:length(data)
    DV_betas{s} = get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', DV_mask);
end


% get RU ROIs from GLM 36
[RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);

% get TU ROIs from GLM 36
[TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);



for c = 1:length(RU_masks)

    disp(RU_masks{c});

    rs = [];
    for s = 1:length(data)
        RU_betas{s} = get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', RU_masks{c});

        [r, p] = corr(RU_betas{s}, DV_betas{s});
        rs(s) = r;
    end

    rs

    rs = atanh(rs);
    [h, p, ci, stats] = ttest(rs)
end
