% compare GLMs using BMS; do CV -- pick ROI for subj s based on contrast for other subjects

clear all;

EXPT = exploration_expt();
[~,~,goodRuns,goodSubjects] = exploration_getSubjectsDirsAndRuns();

data = load_data;

% compare DV vs RT in left M1
%
[masks, region] = get_masks(29, 'DV', true, [], 1);
glms = [29 71 72];

for s = 1:length(data)

    % exclude subj from contrast
    subjs = [1:s-1 s+1:length(data)];
    ccnl_fmri_con(exploration_expt(), 29, {'DV'}, subjs);
    [masks, region] = get_masks(29, 'DV', true, [], 1);

    for i = 1:length(glms)
        % pick top ROI only
        bic = ccnl_bic(EXPT, glms(i), masks{1}, s);
        assert(length(bic) == 1);

        bics(s,i) = bic;
    end
    subj_masks{s} = masks{1};

end

subj_masks'


lme = -0.5 * bics;
[alpha, exp_r, xp, pxp, bor] = bms(lme);

pxp


save('glm_bic_bms_CV.mat');


% revert to original contrast
ccnl_fmri_con(exploration_expt(), 29, ...
    {'DV', 'trial_onset', 'choice_onset', 'feedback_onset'}, ...
     goodSubjects);
