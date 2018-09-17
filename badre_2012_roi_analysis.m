
masks = badre_2012_create_masks(false);

%badre_roi = extract_roi_betas(masks, 'trial_onset');
%save('badre_roi.mat');

load('badre_roi.mat');

for i = 1:length(roi)
    [r, p] = corr(k
end
