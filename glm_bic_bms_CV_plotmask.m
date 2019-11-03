% plot heatmap of masks for different subjects 


load('glm_bic_bms_CV.mat');

conj = [];
EXPT = exploration_expt();

for s = 1:length(subj_masks)
    [mask,V,Y] = ccnl_load_mask(subj_masks{s});
    V.fname = 'glm_bic_bms_CV_heatmap.nii'; % change immediately!!
    if isempty(conj)
        conj = mask;
    else
        conj = conj + mask;
    end
end

spm_write_vol(V, conj);
bspmview(V.fname, fullfile(EXPT.modeldir, 'mean.nii'));
