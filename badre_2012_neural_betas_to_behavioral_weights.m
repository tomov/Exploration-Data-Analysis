% See if RU betas from Badre's 2012 RLPFC ROI correlate with the behavioral weights for RU
%

clear all;

data = load_data;
EXPT = exploration_expt();
glmodel = 21;
regressor = 'RU';

filename = ['badre_2012_b_to_w_glm', num2str(glmodel), '.mat'];
disp(filename);

load results_glme_fig3_nozscore.mat;
w = getEffects(results_VTURU, false);

masks = badre_2012_create_masks(false);

for i = 1:length(masks)
    mask = masks{i};

    clear b;
    for s = 1:length(data)
        b(s) = mean(ccnl_get_beta(EXPT, glmodel, regressor, mask, s));
    end
    all_b{i} = b;

    [r,p] = corr(w(:,2), b');
    disp(mask);
    r
    p

    ps(i,:) = p;
    rs(i,:) = r;
end


p_uncorr = ps;
p_corr = 1 - (1 - ps) .^ length(ps);
r = rs;

save(filename);

table(masks', p_uncorr, p_corr, r);
