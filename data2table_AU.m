function tbl = data2table_AU(roi, data, results)

G = []; N = []; S = []; NL = []; NR = []; GL = []; GR = [];

for s = 1:length(data)
    if size(results.x, 1) == 1
        latents = AU(data(s), results.x);
    else
        latents = AU(data(s), results.x(s,:));
    end
    which = ~data(s).timeout;
    GL = [GL; latents.G(which, 1)];
    GR = [GR; latents.G(which, 2)];
    NL = [NL; latents.N(which, 1)];
    NR = [NR; latents.N(which, 2)];
    G = [G; sum(latents.G(which,:), 2)];
    N = [N; sum(latents.N(which,:), 2)];
    tot = [tot; latents.a * latents.G + latents.b * latents.N];
    n = sum(which);
    S = [S; zeros(n,1) + s];
end

tbl = table(G,N,tot,GL,GR,NL,NR,S);

for roi_idx = 1:length(roi)
    b = [];
    for s = 1:length(data)
        which = ~data(s).timeout;
        b = [b; mean(roi(roi_idx).subj(s).betas(which,:))];

        i = isnan(roi(roi_idx).subj(s).betas(data(s).timeout,:));
        assert(all(i));
    end

    tbl = addvars(tbl, b, 'NewVariableNames', roi(roi_idx).name);
end

