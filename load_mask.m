function [mask, Vmask, Ymask] = load_mask(mask_filename)
% Load a .nii file as a SPM mask
%
% INPUT:
% mask_filename = path to .nii file with mask
%
% OUTPUT:
% mask = binary mask for the non-0 non-nan voxels (from Ymask)
% V = volume struct from spm_vol
% Y = mask from spm_read_vols 
%
Vmask = spm_vol(mask_filename);
Ymask = spm_read_vols(Vmask);
mask = Ymask~=0 & ~isnan(Ymask);
