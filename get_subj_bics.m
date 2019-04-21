function [BICs, logliks] = get_subj_BICs(glme, tbl, exclude)

    y = predict(glme, tbl);

    subjs = unique(tbl.S);
    BICs = [];
    logliks = [];

    p = randomEffects(glme);

    for s = 1:length(subjs)
        S = subjs(s);
        which = tbl.S == S & ~exclude;
        lik = binopdf(tbl.C(which), 1, y(which));
        loglik = sum(log(lik));

        K = length(p) / length(subjs);
        N = sum(which);
        BIC = K * log(N) - 2 * loglik;
        BICs = [BICs; BIC];
        logliks = [logliks; loglik];
    end

end
