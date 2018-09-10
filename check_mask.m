function check_mask( mask , subj)
% check a mask against the mean structural, make sure it looks right
%
% INPUT:
% mask = path of .nii file with mask
%
EXPT = exploration_expt();

if ~exist('subj', 'var')
    struc = fullfile(EXPT.modeldir, 'mean.nii');
else
    struc = fullfile(EXPT.subject(subj).datadir, 'wBrain.nii');
end

P = {struc, fullfile(mask)};
spm_check_registration(char(P));

