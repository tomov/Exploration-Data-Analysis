data = load_data;

fitfile = 'fit_ACU_25nstarts_mixed';
fn = @ACU;
cond = 1; % 1 = RS, 2 = SR, 3 = RR, 4 = SS

load(fitfile, 'results');

figure;

for s = 1:10 %length(data)
    if size(results.x, 1) == 1
        latents = ACU(data(s), results.x);
    else
        latents = ACU(data(s), results.x(s,:));
    end

    G = sum(latents.G(data(s).cond == cond, :), 2);
    N = sum(latents.N(data(s).cond == cond, :), 2);

    subplot(10, 1, s);
    hold on;
    h = plot([G, N]);
    set(h, {'color'}, {[1 0 0]; [0 1 0]});
    for t = 10.5:10:80
        plot([t t], [-1 20], '--', 'color', [0.8 0.8 0.8]);
    end
    legend(h, {'G', 'N'});
    set(gca, 'xtick', []);
    hold off;
end
