[~,~,goodRuns,goodSubjects] = exploration_getSubjectsDirsAndRuns();

% is V in vmPFC?

%{
mni = [ ... %3 30 -21; ... % Daw et al 2006
       6 36 -8 ... % Lim et al 2011
       ];    
       %}
mni = [-3 47 -18];



sphere = 10;
r = sphere / 1.5; % 10 mm sphere

[mask, V, Y] = load_mask('masks/mask.nii');
V.fname = 'masks/badre_vmpfc.nii'; % change immediately!

for i = 1:size(mni, 1)
    cor = mni2cor(mni(i,:), V.mat);

    V.fname = fullfile('masks', sprintf('vmpfc_%d_%d_%d_r=%.1fmm.nii', mni(i,1), mni(i,2), mni(i,3), sphere));

    [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
    sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
    spm_write_vol(V, sphere_mask);

    masks{i} = V.fname;
end









%ccnl_create_mask({'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Rectus_L', 'Rectus_R'}, 'masks/test.nii', 'AAL2');
ccnl_create_mask({'Frontal_Sup_Medial_L', 'Frontal_Med_Orb_L', 'Rectus_L'}, 'masks/vmpfc_l.nii', 'AAL2');


%masks = {'masks/value_association-test_z_FDR_0.01.nii'}
masks = {'masks/vmpfc_l.nii'};

for m = 1:length(masks)

    maskfile = masks{m};

    maskfile

    bs = ccnl_get_beta(exploration_expt, 69, 'V', maskfile, goodSubjects);
    b = nanmean(bs,2);
    b

    [h,p,ci,stat] = ttest(b)
end
