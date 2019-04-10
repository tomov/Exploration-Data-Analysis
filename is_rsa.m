% inter-subject RSA
% from The computational and neural substrates of moral strategies in social decision-making
% Jeroen M. van Baar, Luke J. Chang & Alan G. Sanfey  2019 Nat Comms

clear all;

glmodel = 36;
regressor = 'trial_onset';
EXPT = exploration_expt();

group_mask_filename = fullfile('masks', 'mask.nii');

parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2.nii');
[~, Vparcel, parcellation_vol] = ccnl_load_mask(parcellation_file);
parcellation_vol = round(parcellation_vol);

data = load_data;

parcel_idxs = unique(parcellation_vol(:));


load results_glme_fig3_nozscore.mat;
w = getEffects(results_VTURU, false);

model_RDM = pdist(w, 'euclidean');


for i = 1:length(parcel_idxs)

    parcel_idx = parcel_idxs(i);
    if parcel_idx == 0
        continue;
    end
   
    % ROI
    mask = parcellation_vol == parcel_idx;

    % normalize mask
    group_vol = spm_vol(group_mask_filename);
    group_mask = spm_read_vols(group_vol);

    [x, y, z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
    cor = mni2cor(cor2mni([x y z], Vparcel.mat), group_vol.mat); % voxel coords in AAL2 space --> voxel coords in MNI space --> voxel coords in our space
    ind = sub2ind(size(group_mask), cor(:,1), cor(:,2), cor(:,3)); % voxel coords in our space --> voxel indices
    
    % Reconstruct mask in our space
    %
    Vmask = group_vol;
    Vmask.fname = 'tmp.nii'; % CRUCIAL! o/w overwrite mask.nii
    mask = zeros(size(group_mask));
    mask(ind) = 1; % voxel indices --> binary mask
    
    % Only include voxels that are part of the subject group-level mask
    % i.e. that have non-NaN betas for all subjects
    %
    mask = mask & group_mask;

    % extract betas
    for s = 1:length(data)
        % from ccnl_decode_regressor
        %
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(s)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

        names = SPM.xX.name'; % regressor names
        which_reg = contains(names, regressor);

        % extract betas B
        cdir = pwd;
        cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
        B = spm_data_read(SPM.Vbeta(which_reg), find(mask));
        cd(cdir);

        b(s,:) = nanmean(B, 1);
    end

    all_b{parcel_idx} = b;

    neural_RDM = pdist(b, 'correlation');

    break; % TODO rm
end
