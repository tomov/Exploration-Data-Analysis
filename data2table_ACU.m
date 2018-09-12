function [tbl, lats] = data2table_ACU(roi, data, results)

G = []; N = []; S = []; NL = []; NR = []; GL = []; GR = []; tot = [];

for s = 1:length(data)
    if size(results.x, 1) == 1
        latents = ACU(data(s), results.x);
    else
        latents = ACU(data(s), results.x(s,:));
    end
    lats(s) = latents;

    % TODO figure out if we're excluding timeouts
    GL = [GL; latents.G(:, 1)];
    GR = [GR; latents.G(:, 2)];
    NL = [NL; latents.N(:, 1)];
    NR = [NR; latents.N(:, 2)];
    G = [G; sum(latents.G, 2)];
    N = [N; sum(latents.N, 2)];
    tot = [tot; latents.a * sum(latents.G, 2) + latents.b * sum(latents.N, 2)];
    n = size(latents.G,1);
    S = [S; zeros(n,1) + s];
end

tbl = table(G,N,tot,GL,GR,NL,NR,tot,S);

% notice that we keep the NaN betas for the bad runs
%
for roi_idx = 1:length(roi)
    b = [];
    for s = 1:length(data)
        b = [b; nanmean(roi(roi_idx).subj(s).betas, 2)];
    end

    tbl = addvars(tbl, b, 'NewVariableNames', roi(roi_idx).name);
end

