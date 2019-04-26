% create lateralized copy of neurosynth parcellation

parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2.nii');
[~, Vparcel, parcellation_vol] = ccnl_load_mask(parcellation_file);
parcellation_vol = round(parcellation_vol);

parcel_idxs = unique(parcellation_vol(:));

lateralized_file = fullfile('masks', 'Neurosynth Parcellation_2_lateralized.nii');
lateralized_vol = parcellation_vol;
Vlateral = Vparcel;
Vlateral.fname = lateralized_file; % change immediately!!!!

max_parcel_idx = max(parcel_idxs);
hemi_xs = 1:size(lateralized_vol,1)/2; % one hemisphere
lateralized_vol(hemi_xs,:,:) = lateralized_vol(hemi_xs,:,:) + max_parcel_idx;
lateralized_vol(parcellation_vol == 0) = 0; % clear out 0's

spm_write_vol(Vlateral, lateralized_vol);

ccnl_view_mask({lateralized_file, parcellation_file});
