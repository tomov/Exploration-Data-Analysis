function cross_subj_perf(roi_glmodel, roi_contrast, glmodel, standardize, clusterFWEcorrect, extent, odd_runs)

    % correlate neural BICs with performance 

% TODO dedupe with cross_subject.m

if ~exist('standardize', 'var')
    standardize = false;
end
if ~exist('clusterFWEcorrect', 'var')
    clusterFWEcorrect = true;
end
if ~exist('extent', 'var')
    extent = [];
end
if ~exist('odd_runs', 'var')
    odd_runs = false;  % all runs
end

what = 'sphere';

EXPT = exploration_expt();

data = load_data;

filename = sprintf('cross_subj_perf_bic_roiglm%d_%s_glm%d_%s_standardize=%d_corr=%d_extent=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, what, standardize, clusterFWEcorrect, extent);
disp(filename);


% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent);


perf = getPerfs(data);


for c = 1:length(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        bic(s) = ccnl_bic(EXPT, glmodel, mask, s)
    end
    all_bic{c} = bic;

    [r, p] = corr(perf, bic');

    disp(mask);
    r
    p

    ps(c,:) = p;
    rs(c,:) = r;
end

p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);
r = rs;

save(filename, '-v7.3');

if exist('region', 'var')
    masknames = region';
end
table(masknames', p_uncorr, p_corr, r)
