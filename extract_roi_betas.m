function roi = extract_roi_betas(masks, event)
% extracts betas for a set of masks from the preloaded betas (preload_betas_from_masks for the whole brain)
% much faster than ccnl_get_beta
%
% USAGE:
%   roi = extract_roi_betas(masks, event)

assert(ismember(event, {'trial_onset', 'choice_onset', 'feedback_onset'}));

data = load_data;

whole_brain_mask = 'masks/mask.nii';

[allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();

for s = 1:length(data)

    betas_filename = betapath_from_maskpath(whole_brain_mask, s, event);
    disp(['  betas_filename = ', betas_filename]);
    load(betas_filename, 'b'); % loads b = betas for whole-brain

    for roi_idx = 1:length(masks)
        mask = masks{roi_idx};

        m = load_mask(mask);
        cnt = sum(m(:));

        disp(['    mask ', mask, ': ', num2str(cnt), ' voxels']);
        [~, maskname, ~] = fileparts(mask);

        roi(roi_idx).subj(s).betas = get_activations_submask(m, b);
        roi(roi_idx).mask = mask;
        roi(roi_idx).name = maskname;
    end
end

