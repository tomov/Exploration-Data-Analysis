% extract the betas from the nosmooth GLM 23 for the BG ROIs
% note: betas are NaN for runs / trials that were excluded b/c of motion, etc
% note: betas are included for timeouts!
%

data = load_data;
event = 'trial_onset';
EXPT = exploration_expt();
rsa_glm = 23;

[allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();

rois = {'striatum', 'pallidum', 'caudate', 'putamen'};
clear roi;

for roi_idx = 1:length(rois)
    mask = fullfile('masks', [rois{roi_idx}, '.nii']);
    m = load_mask(mask);
    cnt = sum(m(:));

    for s = 1:length(data)
        runs = find(goodRuns{s});
        b = nan(length(data(s).trial), cnt);

        disp(s);

        for r = 1:length(EXPT.subject(s).functional)
            r = runs(r);
            which = data(s).run == r;
            block = data(s).block(which);
            trial = data(s).trial(which);
            for i = 1:length(block)
                name = [event, '_run_', num2str(r), '_block_', num2str(block(i)), '_trial_', num2str(trial(i))];
                disp(name);
                betas = ccnl_get_beta_nosmooth(EXPT, rsa_glm, name, mask, s);
                b(data(s).run == r & data(s).block == block(i) & data(s).trial == trial(i), :) = betas;
            end
        end

        roi(roi_idx).mask = mask;
        roi(roi_idx).name = rois{roi_idx};
        roi(roi_idx).subj(s).betas = b;
    end
end

save('john_roi.mat', 'roi');

