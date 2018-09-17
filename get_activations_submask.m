function activations = get_activations_submask(my_mask, activations)
% Given a custom mask and a bunch of beta (or t-value) vectors based on the subject group-level
% mask (i.e. get_betas('mask.nii', ...)), return beta (or t-value) vectors corresponding to the custom mask.
% This is in order to avoid having to use ccnl_get_beta for the custom
% mask and to do stuff like spotlight search.
%
% USAGE:
% activations = get_activations_submask(my_mask, activations)
%
% INPUT:
% my_mask = 3D binary mask whose betas (or t-values) we want to extract, in the space of
%           the subject group-level mask
% activations = [nTimePoints x nVoxels] matrix where nVoxels = # voxels in the
%               subject group-level mask (and they correspond to them too)
%
% OUTPUT:
% activations = [nTimePoints x mVoxels] matrix where mVoxels = # voxels in
%         my_mask
%
% EXAMPLES:
% [data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
% whole_brain_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata);
% hippo_mask = load_mask('masks/hippocampus.nii');
% hippo_betas = get_activations_submask(hippo_mask, whole_brain_betas);
%

% Load the subject group-level mask
%
group_mask_filename = fullfile('masks', 'mask.nii'); % the subject group-level mask
[~, group_vol, group_mask] = load_mask(group_mask_filename);

% Make sure custom mask is in same coordinate space as the group-level mask
%
assert(isequal(size(my_mask), size(group_mask)));

% Get the indices of the voxels in the group-level mask; make sure they
% correspond to the voxels in the beta (or t-value) vector
%
group_inds = find(group_mask);
assert(numel(group_inds) == size(activations, 2));

% Get the indices of the voxels in the custom mask
%
my_inds = find(my_mask);

% For each voxel, check if it's part of the custom mask. If it is, then use
% its beta (or t-value).
%
which = ismember(group_inds, my_inds);
activations = activations(:, which);

