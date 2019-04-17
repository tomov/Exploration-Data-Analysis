% inter-subject RSA
% from The computational and neural substrates of moral strategies in social decision-making
% Jeroen M. van Baar, Luke J. Chang & Alan G. Sanfey  2019 Nat Comms

clear all;

printcode;

glmodel = 57;
regressor = 'trial_onset';
EXPT = exploration_expt_nosmooth();
null_iters = 10000;

group_mask_filename = fullfile('masks', 'mask.nii');

parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2.nii');
[~, Vparcel, parcellation_vol] = ccnl_load_mask(parcellation_file);
parcellation_vol = round(parcellation_vol);

data = load_data;

parcel_idxs = unique(parcellation_vol(:));


% compute model RDM
load results_glme_fig3_nozscore.mat;
w = getEffects(results_VTURU, false);
model_RDM = pdist(w, 'cosine');


for i = 1:length(parcel_idxs)

    parcel_idx = parcel_idxs(i);
    if parcel_idx == 0
        continue;
    end

    i
    tic
   
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

    clear b;

    % extract betas
    % TODO ccnl_get_beta_fast
    for s = 1:length(data)
        B = ccnl_get_beta_series(EXPT, glmodel, s, regressor, mask);

        b(s,:) = nanmean(B, 1);
    end

    all_b{parcel_idx} = b;

    % compute neural RDM
    neural_RDM = pdist(b, 'cosine');

    % second-order correlation
    %
    [rho, p] = corr(model_RDM', neural_RDM', 'type', 'Spearman')
    spearman_rho(parcel_idx) = rho;
    spearman_p(parcel_idx) = p;

    toc

    % generate null 
    null_rhos = [];

    disp('Generating null');
    tic

    for iter = 1:null_iters
        % randomly shuffle observations i.e. rows i.e. subjects
        b = b(randperm(size(b,1)), :);

        % recompute neural RDM
        neural_RDM = pdist(b, 'correlation');

        % recompute second-order
        null_rhos(iter) = corr(model_RDM', neural_RDM', 'type', 'Spearman');
    end

    toc

    all_null_rhos{parcel_idx} = null_rhos;

    pval(parcel_idx) = mean(null_rhos >= rho); % what fraction of the null rhos are "better" than rho? = P(we got a value as extreme as rho under the null)

    pval(parcel_idx)

end


save('is_rsa_glmodel=57_trial_onset.mat', '-v7.3');
