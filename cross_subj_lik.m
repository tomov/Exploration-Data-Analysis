function cross_subj_lik(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)

    % correlate neural betas with log likelihood of choices 

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

filename = sprintf('cross_subj_lik_roiglm%d_%s_glm%d_%s_%s_standardize=%d_corr=%d_extent=%d.mat', roi_glmodel, replace(roi_contrast, ' ', '_'), glmodel, regressor, what, standardize, clusterFWEcorrect, extent);
disp(filename);


% get ROIs
[masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent);

% get lik
%

load results_glme_fig3_nozscore.mat;

tbl = data2table(data,0,1); % don't standardize

rmpath('/n/sw/helmod/apps/centos7/Core/spm/12.7487-fasrc01/external/fieldtrip/external/stats/'); % for binopdf; use the MATLAB one, not SPM one

loglik = []; % P(better option)
for s = 1:length(data)
    c = tbl.C(tbl.S == s); % I(a == 1)
    y = predict(results_VTURU, tbl(tbl.S == s,:));  % P(a == 1)
    lik = binopdf(c,1,y);
    loglik = [loglik; sum(log(lik))];
end


for c = 1:length(masks)
    mask = masks{c};
    [~, masknames{c}, ~] = fileparts(mask);

    clear b;
    for s = 1:length(data)
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, ['x', regressor], mask, s));
    end
    all_b{c} = b;

    [r, p] = corr(loglik, b');

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
