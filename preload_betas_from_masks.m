function preload_betas_from_masks(masks, event)
%
% preload betas for all subjects for set of masks 
% masks = cell array of mask paths
%

assert(ismember(event, {'trial_onset', 'choice_onset', 'feedback_onset'}));

data = load_data;
EXPT_nosmooth = exploration_expt_nosmooth();
rsa_glm = 57;

[allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();

for roi_idx = 1:length(masks)
    mask = masks{roi_idx};

    m = load_mask(mask);
    cnt = sum(m(:));

    disp(['mask ', mask, ': ', num2str(cnt), ' voxels']);

    for s = 1:length(data)
        runs = find(goodRuns{s});
        b = nan(length(data(s).trial), cnt);

        disp(['subj ', num2str(s)]);

        for r = 1:length(EXPT_nosmooth.subject(s).functional)
            r = runs(r);
            which = data(s).run == r;
            block = data(s).block(which);
            trial = data(s).trial(which);
            for i = 1:length(block)
                name = [event, '_run_', num2str(r), '_block_', num2str(block(i)), '_trial_', num2str(trial(i))];
                disp(name);
                betas = ccnl_get_beta(EXPT_nosmooth, rsa_glm, name, mask, s);
                b(data(s).run == r & data(s).block == block(i) & data(s).trial == trial(i), :) = betas;
            end
        end

        betas_filename = betapath_from_maskpath(mask, s, event);
        disp(['Saving to ', betas_filename]);
        save(betas_filename, '-v7.3');
    end
end

