load neurosynth_bms_TU_orth=0_standardize=0_mixed=1_intercept=1_method=ridge_getnull=0_zav=0_pa=0_us=0.mat
%load neurosynth_bms_RU_orth=0_standardize=0_mixed=1_intercept=1_method=ridge_getnull=0_zav=0_pa=0_us=0.mat
V = Vparcel;
V.fname = 'masks/neurosynth_bms_RU_BICdiff.nii';



[alpha,exp_r,xp,pxp,bor] = bms(LMEs);

fprintf('BOR = %.6f\n', bor);
fprintf('PXP of original GLM = %.6f\n', pxp(1));
pxp = pxp(2:end);
pxp = pxp';

reg_idx = find(contains(reg_names.Name, 'dec'));

p_uncorr = ps(:,reg_idx);
p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
BIC_orig = BIC(:,1);
BIC_both = BIC(:,2);
p_comp_corr = 1 - (1 - p_comp) .^ numel(p_comp);
T = table(region, p_uncorr, p_corr, pears_rs, pears_ps, BIC_orig, BIC_both, p_comp, p_comp_corr, pxp, p_ax, r_ax);





BIC_diff = BIC_orig - BIC_both;
idx = find(BIC_diff > 0);


mask = parcellation_vol;
mask(:) = 0;

for i = 1:length(idx)
    mask(parcellation_vol == idx(i)) = BIC_diff(idx(i));
end

spm_write_vol(V, mask);

EXPT = exploration_expt();
struc = fullfile(EXPT.modeldir, 'mean.nii');
bspmview(V.fname, struc);


table(idx, BIC_diff(idx))
