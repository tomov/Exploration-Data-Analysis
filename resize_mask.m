function m = resize_mask(inmaskfile, outmaskfile)

    % resize mask into subject native space

    [mask,Vmask] = ccnl_load_mask(inmaskfile);
    [m,V,Y] = ccnl_load_mask('masks/mask.nii'); % native mask
    V.fname = outmaskfile; % change immediately!!!

    [mask_cor1, mask_cor2, mask_cor3] = ind2sub(size(mask),find(mask==1));
    cor = mni2cor(cor2mni([mask_cor1 mask_cor2 mask_cor3], Vmask.mat),V.mat);
    inds = sub2ind(size(Y),cor(:,1),cor(:,2),cor(:,3));


    valid_vox = m;

    m(:) = 0;
    m(inds) = 1;
    m = m & valid_vox; % exclude out-of-brain voxels

    spm_write_vol(V, m);
