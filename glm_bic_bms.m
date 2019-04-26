% compare GLMs using BMS

clear all;

EXPT = exploration_expt();
[~,~,goodRuns,goodSubjects] = exploration_getSubjectsDirsAndRuns();

data = load_data;


[masks, region] = get_masks(47, 'DV', false, 100, 1);

glms = [47 61];

for c = 1:length(masks)
    mask = masks{c};

    mask

    lmes = [];
    bics{c} = [];
    for i = 1:length(glms)
        glmodel = glms(i);
        bic = ccnl_bic(EXPT, glmodel, mask, goodSubjects);
        bics{c} = [bics{c} bic];
    end

    lme = -0.5 * bics{c};
    [alpha, exp_r, xp, pxp, bor] = bms(lme);

    pxp

    pxps(c,:) = pxp;
end

table(region, pxps)