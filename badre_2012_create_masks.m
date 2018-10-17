% create RLPFC masks based on Badre et al. 2012
%
function masks = badre_2012_create_mask(do_create)

% RLPFC from Badre et al. 2012
mni = [36 56 -8; ... % RU @ trial onset
       40 60 -10; ... % explorers; ventral RLPFC
       30 52 -14; ... % same
       35 50 15; ... % Zajkowski et al. 2017 RLPFC
       24 48 20; ... % explorers; dorsal RLPFC
       30 52 16; ... % same
       18 40 22; ... % same
       24 46 20; ... % explorers - nonexplorers
       30 52 -14; ... % RU orth TU, explorers; ventral RLPFC; Fig 4B
       36 56 -10; ... % same
       22 56 26; ... % RU orth TU, explorers; dorsal RLPFC; Fig 4B
       26 52 16; ... % same
       44 42 28; ... % same
       22 54 28; ... % RU orth TU, explore - nonexplore; Fig 4C
       28 48 14; ... % same
       22 46 20];    % same


r = 10 / 1.5; % 10 mm sphere

[mask, V, Y] = load_mask('masks/mask.nii');
V.fname = 'masks/badre_rlpfc.nii'; % change immediately!

for i = 1:size(mni, 1)
    cor = mni2cor(mni(i,:), V.mat);

    V.fname = fullfile('masks', sprintf('badre_rlpfc_%d_%d_%d_r=10mm.nii', mni(i,:)));

    if do_create
        [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
        sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
        spm_write_vol(V, sphere_mask);
    end

    masks{i} = V.fname;
end

%check_mask('masks/badre_rlpfc.nii');
