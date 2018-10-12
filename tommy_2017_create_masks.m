% create RLPFC masks based on Badre et al. 2012
%
function masks = tommy_2017_create_mask(do_create)

% RLPFC from Badre et al. 2012
mni = [-30 16 -8; ... % l insula
       32 22 -8; ... % r insula
       8 16 46];    % dACC


r = 10 / 1.5; % 10 mm sphere

[mask, V, Y] = load_mask('masks/mask.nii');
V.fname = 'masks/tommy.nii'; % change immediately!

regions = {'Insula_L', 'Insula_R', 'dACC_R'};

for i = 1:size(mni, 1)
    cor = mni2cor(mni(i,:), V.mat);

    V.fname = fullfile('masks', sprintf('tommy_%s_%d_%d_%d_r=10mm.nii', regions{i}, mni(i,1), mni(i,2), mni(i,3)));

    if do_create
        [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
        sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
        spm_write_vol(V, sphere_mask);
    end

    masks{i} = V.fname;
end

masks = [masks {'masks/NAC.nii'}];
