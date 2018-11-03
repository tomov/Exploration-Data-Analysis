% create RLPFC masks based on Badre et al. 2012
%
function masks = dlpfc_2012_create_mask(do_create)

% RLPFC from Badre et al. 2012
mni = [40 30 34; ... % TU @ trial onset; neutrally defined DLPFC (from Fig 5)
       38 30 34; ... % from TU @ trial_onset contrast
       30 26 20; ...
       46 14 28];


r = 10 / 1.5; % 10 mm sphere

[mask, V, Y] = load_mask('masks/mask.nii');
V.fname = 'masks/badre_dlpfc.nii'; % change immediately!

for i = 1:size(mni, 1)
    cor = mni2cor(mni(i,:), V.mat);

    V.fname = fullfile('masks', sprintf('badre_dlpfc_%d_%d_%d_r=10mm.nii', mni(i,:)));

    if do_create
        [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
        sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
        spm_write_vol(V, sphere_mask);
    end

    masks{i} = V.fname;
end

%check_mask('masks/badre_dlpfc.nii');
