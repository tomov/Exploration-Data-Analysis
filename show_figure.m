function show_figure(fig)

    switch fig
            
        case 'badre_ROI'
            % show the Badre ROI 
            %
            EXPT = exploration_expt();
            struc = fullfile(EXPT.modeldir,'mean.nii');
            masks = badre_2012_create_masks(false);
            bspmview(masks{1}, struc);


        case 'corr_regs'

            data = load_data;
            tbl = data2table(data,0,1); % don't standardize!

            %data = load_data_sam_2019('data_sam_2019_paper.csv');
            %tbl = data2table_sam_2019(data,0,1); % don't standardize!

            figure;

            regs = {'V', 'RU', 'TU'};

            k = 1;
            for i = 1:3
                for j = i+1:3
                    
                    A = tbl.(regs{i});
                    B = tbl.(regs{j});

                    clear rs;
                    for s = 1:length(data)
                        a = A(tbl.S == s);
                        b = B(tbl.S == s);

                        [r,p] = corr(a, b);
                        rs(s) = r;
                    end

                    k = k+1;

                    m = mean(rs);
                    se = std(rs) / sqrt(length(rs));
                    [h, p, ci, stat] = ttest(atanh(rs));
                    fprintf('%s vs. %s: r = %.2f +- %.2f, t(%d) = %.2f, p = %s\n', regs{i}, regs{j}, m, se, stat.df, stat.tstat, pvalue_to_latex(p));

                end
            end


        case 'fig:simulate'

            figure('pos', [10 10 300 200]);

            fontsize = 12;
            linewidth = 1;
            markersize = 5;

            %
            % compare perf generatively using fitted params
            %

            load simulate.mat

            perf = [perf_V perf_VRU perf_VTU perf_VTURU];
            names = {'softmax (V)', 'UCB (V + RU)', 'Thompson (V/TU)', 'hybrid (V + RU + V/TU)'};
            m = mean(perf,1);
            se = std(perf,[],1) / sqrt(size(perf,1));

            errorbar(1:4, m, se,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            xlabel('Model');
            ylabel('P(better option)');
            xticks(1:4);
            xlim([0 5]);
            ylim([0.73 0.83]);
            xticklabels(names);
            xtickangle(30);

            print('images/simulate', '-dpdf');



            [p, anovatab, stats] = anova1(perf);
            [c,m,h,nms] = multcompare(stats);

            %comp_names = {'V vs. V + RU', 'V vs. V/TU', 'V + RU vs. V + RU + V/TU', 'V/TU vs. V + RU + V/TU'};
            comp_names = {'softmax vs. UCB', 'softmax vs. Thompson', 'UCB vs. hybrid', 'Thompson vs. hybrid'};
            comps = [1 2; 1 3; 2 4; 3 4];

            fprintf('significant difference between models ($F(%d,%d) = %.2f, %s$, one-way ANOVA)\n', anovatab{2,3}, anovatab{3,3}, anovatab{2,5}, pvalue_to_latex(anovatab{2,6}));

            for i = 1:4
                j = find(c(:,1) == comps(i,1) & c(:,2) == comps(i,2));
                p_string = pvalue_to_latex(c(j, end));

                fprintf('%s ($%s$, pairwise multiple comparison test)\n', comp_names{i}, p_string);
            end


        case 'fig:simulate_more'

            %
            % compare perf generatively using different params
            %

            load results_glme_fig3_nozscore.mat;

            %load simulate_more.mat 
            load simulate_more_nws=16_niters=10.mat;

            subj_w1s = logical(size(w1s));
            subj_w2s = logical(size(w2s));
            subj_w3s = logical(size(w3s));

            w = getEffects(results_VTURU, true);
            for s = 1:size(w,1)
                [~, i1(s)] = min(abs(w1s - w(s,1)));
                [~, i2(s)] = min(abs(w2s - w(s,2)));
                [~, i3(s)] = min(abs(w3s - w(s,3)));
            end


            figure('pos', [10 10 700 600]);

            xtix = {};
            for i = 1:nws
                xtix{i} = sprintf('%.2f\n', w2s(i));
                ytix{i} = sprintf('%.2f\n', w3s(i));
            end

            col = colormap;
            for i = 1:16
                subplot(4,4,i);

                img = squeeze(m(i,:,:));
                imagesc(img);
                set(gca, 'xtick', []);
                set(gca, 'ytick', []);
                xtickangle(90);
                if ismember(i, [13 14 15 16])
                    xlabel('w_2');
                    xticklabels(xtix(1:2:end));
                    xticks(1:nws);
                end
                if ismember(i, [1 5 9 13])
                    ylabel('w_3');
                    yticklabels(ytix(1:2:end));
                    yticks(1:nws);
                end
                title(sprintf('w_1 = %.2f', w1s(i)));
                colorbar;

                for s = 1:size(w,1)
                    if i == i1(s)
                        % circle where the peeps are
                        hold on;
                        draw_circle(i2(s), i3(s), 1, 'red');
                        hold off;
                    end
                end

                set(gca, 'YDir', 'normal');
            end

            save shit.mat

            print('images/simulate_more', '-dpdf');

        case 'tab:compare'

            data = load_data;
            load results_glme_fig3_nozscore.mat;
            w = getEffects(results_VTURU, false);

            comp{1} = compare(results_V, results_VRU);
            comp{2} = compare(results_V, results_VTU);
            comp{3} = compare(results_VRU, results_VTURU);
            comp{4} = compare(results_VTU, results_VTURU);

            fprintf('\n\n\n');
            fprintf('\\textbf{Model} & \\textbf{Regressors} & \\textbf{AIC} & \\textbf{BIC}  & \\textbf{LL}  & \\textbf{Deviance}  \\\\ \\hline\n');
            fprintf('\n\n\n');

            model_regs = {'V', 'V + RU', 'V/TU', 'V + RU + V/TU'};
            model_names = {'Softmax', 'UCB', 'Thompson sampling', 'UCB/Thompson hybrid'};
            BIC = [comp{1}.BIC(1) comp{3}.BIC(1) comp{4}.BIC(1) comp{4}.BIC(2)];
            AIC = [comp{1}.AIC(1) comp{3}.AIC(1) comp{4}.AIC(1) comp{4}.AIC(2)];
            LogLik = [comp{1}.LogLik(1) comp{3}.LogLik(1) comp{4}.LogLik(1) comp{4}.LogLik(2)];
            Dev  = [results_V.ModelCriterion.Deviance results_VRU.ModelCriterion.Deviance  results_VTU.ModelCriterion.Deviance  results_VTURU.ModelCriterion.Deviance];

            for i = 1:length(BIC)
                fprintf(' %s &  %s & %.2f & %.2f & %.2f & %.2f \\\\ \\hline \n', model_names{i}, model_regs{i}, BIC(i), AIC(i), LogLik(i), Dev(i));
            end


            fprintf('\n\n\n');
            fprintf('\\textbf{Comparison} &  \\textbf{LR-stat} & \\textbf{p-value} \\\\ \\hline \n');
            fprintf('\n\n\n');

            comp_names = {'V vs. V + RU', 'V vs. V/TU', 'V + RU vs. V + RU + V/TU', 'V/TU vs. V + RU + V/TU'};

            for i = 1:length(comp)
                p_string = pvalue_to_latex(comp{i}.pValue);

                fprintf(' %s & %.2f  & %s \\\\ \\hline \n', comp_names{i}, comp{i}.LRStat(2), p_string);
            end


        case 'fig:perf'

            % do subjects with greater w's also perform better?
            %

            data = load_data;
            load results_glme_fig3_nozscore.mat;
            w = getEffects(results_VTURU, false);

            perf = [];
            for s = 1:length(data)
                which = ~data(s).timeout;
                better = (data(s).mu2(which) > data(s).mu1(which)) + 1; 
                C = double(data(s).choice(which) == better); % human choices
                %perf = [perf; mean(data(s).reward)];
                perf = [perf; mean(C)];
            end

            [r, p] = corr(w(:,1), perf);
            fprintf('performance and w_1 ($r = %.2f, p = %.4f$, Pearson correlation across %d subjects)\n', r, p, length(perf));

            [r, p] = corr(w(:,2), perf);
            fprintf('performance and w_2 ($r = %.2f, p = %.3f$, Pearson correlation across %d subjects)\n', r, p, length(perf));

            [r, p] = corr(w(:,3), perf);
            fprintf('performance and w_3 ($r = %.2f, p = %.3f$, Pearson correlation across %d subjects)\n', r, p, length(perf));

            figure('pos', [10 10 700 200]);

            titls = {'softmax', 'UCB', 'Thompson sampling'};
            for i = 1:3
                subplot(1,3,i);
                scatter(w(:,i), perf);
                lsline;
                title(titls{i});
                xlabel(sprintf('w_%d', i));
                ylabel('P(better option)');
            end

            print('images/perf', '-dpdf');


        case 'fig:vifs'

            figure('pos', [10 10 600 600]);

            fontsize = 12;
            linewidth = 3;
            markersize = 3;
            %{
            EXPT = exploration_expt();
            [vifs, names] = ccnl_vifs(EXPT, 36);

            save vifs.mat;
            %}

            load vifs.mat;

            regs = {'RU', 'TU', 'V', 'V/TU'};
            suffix = {'xRU', 'xTU', 'xV^', 'xVTU'};

            grid = zeros(17, 9);
            grid(1:8, 1:4) = reshape(1:32, [4 8])';
            grid(1:8, 6:9) = reshape((1:32) + 32 * 1, [4 8])';
            grid(10:17, 1:4) = reshape((1:32) + 32 * 2, [4 8])';
            grid(10:17, 6:9) = reshape((1:32) + 32 * 3, [4 8])';

	        %[ha, pos] = tight_subplot(17,9,[.01 .03],[.1 .01],[.01 .01])

            for i = 1:length(regs)

                all = [];
                for s = 1:length(vifs)
                    reg_vifs = vifs{s}(contains(names{s}, suffix{i}));
                    all = [all reg_vifs];

                    idx = s + 32 * (i - 1);
                    idx = find(grid == idx);
                    [r, c] = ind2sub(size(grid), idx);

                    %axes(ha(c + (r - 1) * 9)); 
                    ax = subplot(17, 9, c + (r - 1) * 9);
                    ax.Position(2) = ax.Position(2) - 0.0001;

                    hold on;
                    for j  = 1:length(reg_vifs)
                        if reg_vifs(j) > 10
                            plot(j, reg_vifs(j), 'ko', 'MarkerFaceColor', [1 0.5 0], 'markersize', markersize);
                        else
                            plot(j, reg_vifs(j), 'ko', 'MarkerFaceColor', [0.5 1 0], 'markersize', markersize);
                        end
                    end
                    plot([0 length(reg_vifs)], [1 1], 'k');
                    %plot([0 length(reg_vifs)], [2 2], 'b--');
                    %plot([0 length(reg_vifs)], [4 4], 'r--');
                    %plot([0 length(reg_vifs)], [8 8], 'r-');
                    plot([0 length(reg_vifs)], [10 10], 'r-');
                    set(gca, 'xtick', []);
                    set(gca, 'fontsize', 6);
                    hold off;

                    if s == 2
                        title(['     ', regs{i}, ' VIFs for each subject'], 'fontsize', fontsize);
                    end
                    if s == 30
                      %  xlabel('run');
                    end
                end
                fprintf('# of VIFs > 10 for %s: %d out of %d (%.3f %%)\n', regs{i}, sum(all > 10), numel(all), sum(all > 10) / numel(all) * 100);
            end


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.50, 0.95, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.50, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.50, 0.50, 'D', 'FontSize', 25, 'FontWeight', 'bold');

            h = gcf;
            set(h,'PaperOrientation','landscape');
            print('images/vifs', '-dpdf');



        case 'fig:recovery'

            figure('pos', [10 10 520 350]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;


            %load recovery.mat % expo [0, infty]
            %load recovery_mvnrnd.mat % mvnrnd <---- from preprint
            load recovery_gen_mvnrnd.mat % mvnrnd, generative (i.e. generate choices too; don't use subject choices)

            for i = 1:3
                [r,p] = corr(w_orig(:,i), w_rec(:,i));
                fprintf('recovery w_%d: r = %.4f, p = %e\n', i, r, p);

                subplot(2,3,i);
                scatter(w_orig(:,i), w_rec(:,i));
                title(sprintf('w_%d', i));
                xlabel('simulated');
                ylabel('fitted');
                xlim([-5 5]);
                ylim([-5 5]);
                lsline;
            end

            k = 0;
            for i = 1:3
                for j = i+1:3
                    k = k + 1;

                    [r,p] = corr(w_rec(:,i), w_rec(:,j));
                    fprintf('recovery w_%d vs w_%d: r = %.4f, p = %f\n', i, j, r, p);

                    subplot(2,3,k + 3);
                    scatter(w_rec(:,i), w_rec(:,j));
                    title(sprintf('w_%d vs. w_%d', i, j));
                    xlabel(sprintf('fitted w_%d', i));
                    ylabel(sprintf('fitted w_%d', j));
                    xlim([-5 5]);
                    ylim([-5 5]);
                    lsline;
                end
            end

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.50, 'B', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/recovery', '-dpdf');


        case 'fig:learning'

            rng default;

            figure('pos', [10 10 420 450]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;


            conds = {'RS', 'SR', 'RR', 'SS'};

            captions = {'human', 'model'};

            %  learning curves
            %
            data = load_data;
            tbl = data2table(data,0,1); % don't standardize! TODO never standardize anywhere
            load results_glme_fig3_nozscore.mat;
            y = predict(results_VTURU, tbl); % hybrid model predictions

            for human_or_model = 1:2
                for cond = 1:4
                    latents = kalman_filter(data(1));
                    v = linspace(min(latents(1).m(:)),max(latents(1).m(:)),8)';

                    %b = [];
                    for s = 1:length(data)
                        latents = kalman_filter(data(s));

                        which = ~data(s).timeout;

                        if human_or_model == 1
                            better = (data(s).mu2(which) > data(s).mu1(which)) + 1; 
                            C = double(data(s).choice(which) == better); % human choices
                        else
                            better = data(s).mu1(which) > data(s).mu2(which); 
                            C = double(binornd(1, y(tbl.S == s)) == better); % model choices
                        end

                       
                        for t = 1:max(data(s).trial)
                            ix = data(s).trial(which) == t & data(s).cond(which) == cond;
                            if ~any(ix)
                                pc(s,t,cond) = nan;
                            else
                                pc(s,t,cond) = nanmean(C(ix));
                            end
                        end
                    end
                end

                subplot(2,1, human_or_model);

                [se,mu] = wse(pc);
                x = 1:max(data(s).trial);
                errorbar(x,mu(:,1),se(:,1),'-ok','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor','k'); hold on
                errorbar(x,mu(:,2),se(:,2),'-o','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);
                errorbar(x,mu(:,3),se(:,3),'-o','LineWidth',linewidth,'MarkerSize',markersize);
                errorbar(x,mu(:,4),se(:,4),'-o','LineWidth',linewidth,'MarkerSize',markersize);
                legend(conds,'FontSize',fontsize,'Location','SouthEast');
                %set(gca,'FontSize',fontsize,'XLim',[min(v) max(v)],'YLim',[0 1]);
                ylabel('P(better option)','FontSize',fontsize);
                xlabel('trial','FontSize',fontsize);
                title(captions{human_or_model}, 'fontsize', fontsize);

            end

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.50, 'B', 'FontSize', 25, 'FontWeight', 'bold');


            tbl = data2table(data,0,1);
            y = predict(results_VTURU, tbl);
            %mu = [tbl.mu1 tbl.mu2];
            r = [tbl.r1 tbl.r2];
            c = (y < 0.5) + 1;

            for cond = 1:4
                r_cond = r(double(tbl.cond) == cond,:);
                c_cond = c(double(tbl.cond) == cond,:);
                ix = sub2ind(size(r_cond), 1:size(r_cond,1), c_cond'); % chosen r's

                fprintf('%s -> avg model reward = %.3f\n', conds{cond}, mean(r_cond(ix)));
            end

            print('images/learning', '-dpdf');





        case 'fig:psycho'

            figure('pos', [10 10 420 650]);

            fontsize = 12;
            linewidth = 1;
            markersize = 5;


            conds = {[1 2], [3 4]};
            legends = {{'RS', 'SR'}, {'RR', 'SS'}};

            captions = {'human', 'model'};

            %  psychometric curves
            %

            for human_or_model = 1:2
                for curve = 1:2

                    clear pc;
                    clear pc_fit;
                
                    data = load_data;
                    tbl = data2table(data,0,1);

                    load results_glme_fig3_nozscore.mat;
                    y = predict(results_VTURU, tbl); % hybrid model predictions

                    load results_glme_fig4_nozscore.mat
                    yh_probit = predict(results, tbl); % probit for humans

                    % same as Figure2 but for model
                    Cm = binornd(1, y); % model choices
                    tbl = [tbl table(Cm)];
                    if ~exist('results_glme_learning_nozscore.mat', 'file')
                        formula = 'Cm ~ -1 + cond + cond:V + (-1 + cond + cond:V|S)';
                        % note: Laplace method not used here because it doesn't seem to complete
                        results = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','DummyVarCoding','Full');
                        save results_glme_learning_nozscore results
                    else
                        load results_glme_learning_nozscore.mat;
                    end
                    ym_probit = predict(results, tbl); % probit for hybrid model


                    latents = kalman_filter(data(1));
                    v = linspace(min(latents(1).m(:)),max(latents(1).m(:)),8)';
                    %v = linspace(-25, 25, 11)';

                    %b = [];
                    for s = 1:length(data)
                        latents = kalman_filter(data(s));

                        which = ~data(s).timeout;
                        V = latents.m(which,1) - latents.m(which,2);

                        if human_or_model == 1
                            C = double(data(s).choice(which)==1); % human choices
                            C_fit = yh_probit(tbl.S == s);
                        else
                            C = binornd(1, y(tbl.S == s)); % model choices
                            C_fit = ym_probit(tbl.S == s);
                        end

                        
                        for j = 1:length(v)-1
                            ix = V>v(j) & V<v(j+1) & data(s).cond(which) == conds{curve}(1);
                            if ~any(ix)
                                pc(s,j,1) = nan;
                                pc_fit(s,j,1) = nan;
                            else
                                pc(s,j,1) = nanmean(C(ix));
                                pc_fit(s,j,1) = nanmean(C_fit(ix));
                            end
                            
                            ix = V>v(j) & V<v(j+1) & data(s).cond(which) == conds{curve}(2);
                            if ~any(ix)
                                pc(s,j,2) = nan;
                                pc_fit(s,j,2) = nan;
                            else
                                pc(s,j,2) = nanmean(C(ix));
                                pc_fit(s,j,2) = nanmean(C_fit(ix));
                            end
                        end
                    end

                    subplot(4,2, 1 + human_or_model - 1 + 2 * (curve - 1));

                    [se,mu] = wse(pc);
                    [se_fit,mu_fit] = wse(pc_fit);

                    x = v(1:end-1) + diff(v)/2;

                    hold on;
                    p1 = plot(x,mu_fit(:,1),'LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor','k', 'Color', [0 0 0]);
                    errorbar(x,mu(:,1),se(:,1),'ok','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor','k', 'Color', [0 0 0]);

                    p2 = plot(x,mu_fit(:,2),'LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor','k', 'Color',[0.5 0.5 0.5]);
                    errorbar(x,mu(:,2),se(:,2),'o','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);

                    legend([p1 p2], legends{curve},'FontSize',fontsize,'Location','SouthEast');
                    set(gca,'FontSize',fontsize,'XLim',[min(v) max(v)],'YLim',[0 1]);
                    ylabel('Choice probability','FontSize',fontsize);
                    xlabel('Expected value difference','FontSize',fontsize);
                    if curve == 1
                        title(captions{human_or_model});
                    end
                end
            end


            % same as Figure 2 TODO dedupe
            load results_glme_fig4_nozscore.mat
            
            % hypothesis tests
            H = [1 -1 0 0 0 0 0 0; 0 0 1 -1 0 0 0 0; 0 0 0 0 1 -1 0 0; 0 0 0 0 0 0 1 -1];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. SR: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, RR vs. SS: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['slope, RS vs. SR: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['slope, RR vs. SS: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. 0: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, SR vs. 0: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['intercept, RR vs. 0: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['intercept, SS vs. 0: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            % elot results
            [beta,~,stats] = fixedEffects(results);
            subplot(4,2,5);
            errorbar(beta(1:4),stats.SE(1:4),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(1:4),(stats.Upper(1:4) - stats.Lower(1:4))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5],'YLim', [-1 1]);
            ylabel('Intercept','FontSize',fontsize);
            ylim([-0.5 0.5]);

            subplot(4,2,7);
            errorbar(beta(5:8),stats.SE(5:8),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(5:8),(stats.Upper(5:8) - stats.Lower(5:8))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [0 3]);
            ylabel('Slope','FontSize',fontsize);
            ylim([0.075 0.2]);




            % same as Figure 2 but for model TODO dedupe
            load results_glme_learning.mat;
            
            % hypothesis tests
            H = [1 -1 0 0 0 0 0 0; 0 0 1 -1 0 0 0 0; 0 0 0 0 1 -1 0 0; 0 0 0 0 0 0 1 -1];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. SR: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, RR vs. SS: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['slope, RS vs. SR: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['slope, RR vs. SS: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. 0: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, SR vs. 0: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['intercept, RR vs. 0: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['intercept, SS vs. 0: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            % elot results
            [beta,~,stats] = fixedEffects(results);
            subplot(4,2,6);
            errorbar(beta(1:4),stats.SE(1:4),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(1:4),(stats.Upper(1:4) - stats.Lower(1:4))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [-0.3 0.5]);
            ylabel('Intercept','FontSize',fontsize); 
            %ylim([-0.5 0.5]);

            subplot(4,2,8);
            errorbar(beta(5:8),stats.SE(5:8),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(5:8),(stats.Upper(5:8) - stats.Lower(5:8))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [1.2 2.3]);
            ylabel('Slope','FontSize',fontsize);
            %ylim([1 2.5]);


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.03, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.03, 0.50, 'B', 'FontSize', 25, 'FontWeight', 'bold');


            print('images/psycho', '-dpdf');




        case 'fig:task' % Figure1
            figure('pos', [10 10 520 450]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;

            h = subplot(2,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 0.6;
            pos(2) = pos(2) * 0.93;
            pos(3) = pos(3) * 1.2;
            pos(4) = pos(4) * 1.2;
            subplot(2,2, 1, 'position', pos);
            imshow('images/figure1Aleft.png'); %, 'InitialMagnification', 'fit');  


            %% ------------------- Sam_Figure1
            subplot(2,2,2);

            x = linspace(-10,10,100);
            plot(x,normpdf(x,2,4),'-k','LineWidth',linewidth);
            hold on;
            x = -1; y = normpdf(x,x,4);
            plot([x x],[0 y],'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5])
            text(2-0.4,y+0.01,'R','FontSize',fontsize);
            text(-1-0.4,y+0.01,'S','FontSize',fontsize);
            set(gca,'FontSize',fontsize,'YLim',[0 0.12]);
            ylabel('Probability density','FontSize',fontsize);
            xlabel('Reward','FontSize',fontsize);

            %% ---------------- Sam_Figure2

            v = linspace(-30,30,8);
            x = v(1:end-1) + diff(v)/2;
           
            data = load_data;
            tbl = data2table(data,0,1);
            load results_glme_fig3_nozscore.mat;
            CC = predict(results_VTURU,tbl);
            
            for s = 1:length(data)
                latents = kalman_filter(data(s));
                V = latents.m(~data(s).timeout,1) - latents.m(~data(s).timeout,2);
                C = double(data(s).choice(~data(s).timeout)==1);
                c = CC(tbl.S==s);
                
                for j = 1:length(v)-1
                    for cond = 1:4
                        ix = V>v(j) & V<v(j+1) & data(s).cond(~data(s).timeout)==cond;
                        if ~any(ix)
                            pc0(s,j,cond) = nan;
                            pc1(s,j,cond) = nan;
                        else
                            pc0(s,j,cond) = nanmean(C(ix));
                            pc1(s,j,cond) = nanmean(c(ix));
                        end
                    end
                end
            end
            
            mu = linspace(-3,3,100);
            d = 1;
            p(1,:) = 1-normcdf(0,mu+d,1);
            p(2,:) = 1-normcdf(0,mu,1+d);
            
            T = {{'Relative uncertainty: intercept shift'}, {'Total uncertainty: slope shift'}};
            for i = 1:2
                subplot(2,2,i+2);
                plot(mu,p(i,:),'-k','LineWidth',linewidth); hold on;
                plot(mu,1-normcdf(0,mu,1),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                if i==1
                    legend({'RS' 'SR'},'FontSize',fontsize,'Location','southeast');
                else
                    legend({'RR' 'SS'},'FontSize',fontsize,'Location','southeast');
                end
                set(gca,'FontSize',fontsize,'XLim',[min(mu) max(mu)],'YLim',[-0.05 1.05]);
                ylabel('P(choose 1)','FontSize',fontsize);
                xlabel({'Expected value difference,', '\mu(1) - \mu(2)'},'FontSize',fontsize);
                title(T{i},'FontSize',fontsize','FontWeight','Bold');
            end
           

            %{
            [se,mu] = wse(pc0);
            mu_c = squeeze(nanmean(pc1));
            
            for i = 1:2
                subplot(4,2,i+4);
                
                if i==1
                    p1 = errorbar(x/10,mu(:,1),se(:,1),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,2),se(:,2),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,1),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,2),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RS' 'SR'},'FontSize',fontsize,'Location','East');
                else
                    p1 = errorbar(x/10,mu(:,3),se(:,3),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,4),se(:,4),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,3),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,4),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RR' 'SS'},'FontSize',fontsize,'Location','East');
                end
                
                set(gca,'FontSize',fontsize,'XLim',[-3 3],'YLim',[-0.05 1.05]);
                ylabel('P(choose left)','FontSize',fontsize);
                xlabel({'Expected value difference', '(V = Q(left) - Q(right))'},'FontSize',fontsize);
            end
            %}
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.50, 0.95, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.50, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.50, 0.50, 'D', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/task', '-dpdf');




        case 'fig:behav' % Figure2
            figure('pos', [10 10 700 200]);

            % probit regression results

            fontsize = 12;
            linewidth = 3;
            markersize = 6;

            %% -------------------- Sam Figure 4

            % Probit analysis of conditions
            data = load_data;
            
            % fit generalized linear mixed effects model
            if ~exist('results_glme_fig4_nozscore.mat', 'file')
                tbl = data2table(data,0,1);
                formula = 'C ~ -1 + cond + cond:V + (-1 + cond + cond:V|S)';
                % note: Laplace method not used here because it doesn't seem to complete
                results = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','DummyVarCoding','Full');
                save results_glme_fig4_nozscore results
            else
                load results_glme_fig4_nozscore.mat
            end
            
            % hypothesis tests
            H = [1 -1 0 0 0 0 0 0; 0 0 1 -1 0 0 0 0; 0 0 0 0 1 -1 0 0; 0 0 0 0 0 0 1 -1];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. SR: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, RR vs. SS: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['slope, RS vs. SR: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['slope, RR vs. SS: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. 0: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, SR vs. 0: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['intercept, RR vs. 0: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['intercept, SS vs. 0: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);
            
            % elot results
            [beta,~,stats] = fixedEffects(results);
            subplot(1,3,1);
            errorbar(beta(1:4),stats.SE(1:4),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(1:4),(stats.Upper(1:4) - stats.Lower(1:4))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5],'YLim', [-1 1]);
            ylabel('Intercept','FontSize',fontsize);
            ylim([-0.5 0.5]);

            subplot(1,3,2);
            errorbar(beta(5:8),stats.SE(5:8),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta(5:8),(stats.Upper(5:8) - stats.Lower(5:8))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [0 3]);
            ylabel('Slope','FontSize',fontsize);
            ylim([0.075 0.2]);

            % TODO nozscore
            %{
            subplot(1,3,3);

            % Probit analysis of computational variables
            load results_glme_fig3
            results = results_VTURU;
                        
            % plot results
            [beta,~,stats] = fixedEffects(results);
            errorbar(beta([3 1 2]),stats.SE([3 1 2]),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta([3 1 2]),(stats.Upper([3 1 2]) - stats.Lower([3 1 2]))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick',[1 2 3],'XTickLabel',{'$V$' '$RU$' '$V/TU$'},'XLim',[0.5 3.5], 'Ylim', [0 3]);
            ylabel('Regression coefficient','FontSize',fontsize);
            %}

            load results_glme_fig3_nozscore; % same as TrustRegion2D
            results = results_VTURU;
            [beta,~,stats] = fixedEffects(results);
            w = beta([3 1 2]);
            se = stats.SE([3 1 2]);
            t = stats.tStat([3 1 2]);
            df = stats.DF([3 1 2]);
            p = stats.pValue([3 1 2]);
           
            for i = 1:3
                fprintf('w_%d = %.3f \\pm %.3f, t(%d) = %.2f, p = %s; \n', i, w(i), se(i), df(i), t(i), pvalue_to_latex(p(i)));
            end
            
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.07, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.35, 0.95, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            %text(0.64, 0.95, 'C', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/behav', '-dpdf');





        case 'tab:RU'

            %ccnl_results_table('AAL2', 'peak', exploration_expt, 45, 'RU', 0.001, '-', 0.05, 20, 3, true);
            %ccnl_results_table('AnatomyToolbox', 'peak', exploration_expt, 45, 'RU', 0.001, '-', 0.05, 20, 3, true);
            ccnl_results_table('HarvardOxford-maxprob-thr0', 'peak', exploration_expt, 45, 'RU', 0.001, '-', 0.05, 20, 3, true);

        case 'tab:TU'

            ccnl_results_table('AAL2', 'peak', exploration_expt, 45, 'TU', 0.001, '+/-', 0.05, 20, 3, true);
            %ccnl_results_table('AnatomyToolbox', 'peak', exploration_expt, 45, 'TU', 0.001, '+/-', 0.05, 20, 3, true);
            %ccnl_results_table('HarvardOxford-maxprob-thr0', 'peak', exploration_expt, 45, 'TU', 0.001, '+/-', 0.05, 20, 3, true);

        case 'tab:DV'

            %ccnl_results_table('AAL2', 'peak', exploration_expt, 29, 'DV', 0.001, '-', 0.05, 20, 3, true);
            ccnl_results_table('AnatomyToolbox', 'peak', exploration_expt, 29, 'DV', 0.001, '-', 0.05, 20, 3, true);
            %ccnl_results_table('HarvardOxford-maxprob-thr0', 'peak', exploration_expt, 45, 'TU', 0.001, '+/-', 0.05, 20, 3, true);


        case 'fig:RU'

            % RU contrast  from GLM 45
            %
            figure('pos', [100 100 650 180]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(1,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.60;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(2,1, 1, 'position', pos);

            %PICpng = imread('images/badre_RLPFC.png');   %  <-- to cross-check ; bspmview('masks/badre_rlpfc_36_56_-8_r=10.0mm.nii', '../glmOutput/mean.nii')
            %PICpng = imread('images/RU-trial.png');
            %PICpng = imread('images/RU-trial_100.png'); % extent >= 100
            PICpng = imread('images/RU_uncorr.png'); % extent >= 100

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  

            % Badre ROI
            centx = x * 0.795;
            centy = y * 0.29;

            % our ROI
            %centx = x * 0.79;
            %centy = y * 0.30;

            r = 25;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            title('Relative uncertainty (uncorr.)', 'FontSize', fontsize);


            subplot(1,4,3);

            RU_roi_idx = 1;

            load('main_effect_roiglm-1_badre_glm45_RU_corr=0_extent=100_Num=1_s=10.0.mat');
            beta(1) = m(RU_roi_idx);
            ci(1) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(1) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{1} = bs{RU_roi_idx};
            pp(1) = p_uncorr(RU_roi_idx);

            load('main_effect_roiglm-1_badre_glm45_TU_corr=0_extent=100_Num=1_s=10.0.mat');
            beta(2) = m(RU_roi_idx);
            ci(2) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(2) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{2} = bs{RU_roi_idx};
            pp(2) = p_uncorr(RU_roi_idx);


            % t-tests
            %
            [h, p, ci, stats] = ttest(betas{1});
            fprintf('t-test RU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            [h, p, ci, stats] = ttest(betas{2});
            fprintf('t-test TU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            [h, p, ci, stats] = ttest(betas{1}, betas{2});
            fprintf('paired t-test RU - TU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            %for s = 1:length(betas{1})
            %    plot([1 2], [betas{1}(s) betas{2}(s)], '-o', 'MarkerSize', 2,'MarkerFaceColor','k', 'Color', 'k');
            %end

            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
            end
            y = max(beta + err + 0.06);
            line([1 2], [y y], 'color', 'black');
            text(1.45, y + 0.02, significance(p));
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.05 0.25]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title({'Main effect', 'RLPFC (R) [36 56 -8]'}, 'FontSize', axisfontsize);



            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.08, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/RU', '-dpdf');










        case 'fig:TU'

            % TU contrast  from GLM 45
            %
            figure('pos', [100 100 650 180]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(1,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.60;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(2,1, 1, 'position', pos);

            %PICpng = imread('images/badre_DLPFC.png');   %  <-- to cross-check  ; bspmview('masks/badre_dlpfc_38_30_34_r=10.0mm.nii', '../glmOutput/mean.nii')
            %PICpng = imread('images/TU-trial.png');
            %PICpng = imread('images/TU-trial_100.png'); % extent >= 100
            PICpng = imread('images/TU_uncorr.png'); % extent >= 100

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  

            % Badre ROI
            centx = x * 0.76;
            centy = y * 0.15;

            % our ROI
            %centx = x * 0.79;
            %centy = y * 0.30;

            r = 25;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            title('Total uncertainty (uncorr.)', 'FontSize', fontsize);


            subplot(1,4,3);

            TU_roi_idx = 2;

            load('main_effect_roiglm-1_dlpfc_glm45_RU_corr=0_extent=100_Num=1_s=10.0.mat');
            beta(1) = m(TU_roi_idx);
            ci(1) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
            err(1) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
            betas{1} = bs{TU_roi_idx};
            pp(1) = p_uncorr(TU_roi_idx);

            load('main_effect_roiglm-1_dlpfc_glm45_TU_corr=0_extent=100_Num=1_s=10.0.mat');
            beta(2) = m(TU_roi_idx);
            ci(2) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
            err(2) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
            betas{2} = bs{TU_roi_idx};
            pp(2) = p_uncorr(TU_roi_idx);


            % t-tests
            %
            [h, p, ci, stats] = ttest(betas{1});
            fprintf('t-test RU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            [h, p, ci, stats] = ttest(betas{2});
            fprintf('t-test TU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            [h, p, ci, stats] = ttest(betas{2}, betas{1});
            fprintf('paired t-test TU - RU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            %for s = 1:length(betas{1})
            %    plot([1 2], [betas{1}(s) betas{2}(s)], '-o', 'MarkerSize', 2,'MarkerFaceColor','k', 'Color', 'k');
            %end

            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
            end
            y = max(beta + err + 0.06);
            line([1 2], [y y], 'color', 'black');
            text(1.45, y + 0.02, significance(p));
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.05 0.25]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title({'Main effect', 'DLPFC (R) [38 30 34]'}, 'FontSize', axisfontsize);



            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.08, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/TU', '-dpdf');



        case 'tab:augmented'

            fprintf('\\textbf{Model} & \\textbf{Regressors} & \\textbf{AIC} & \\textbf{BIC}  & \\textbf{LL}  & \\textbf{Deviance}  \\\\ \\hline\n');

            %load results_glme_fig3_nozscore.mat; TODO make them the same, fix everywhere
            % note that we use results_orig from the decoder, which EXCLUDES bad_runs
            % but for behavior we use all runs for more power => we're legit
            % just clarify it
            %

            % behavioral GLM
            load('univariate_decoder_refactored_roiglm-1_badre_glm45_RU_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_orig');
            fprintf('UCB/Thompson hybrid       &  V + RU + V/TU  &  %.2f  &  %.2f  &  %.2f  & %.2f    \\\\  \\hline \\hline\n', ...
                results_orig.ModelCriterion.AIC, results_orig.ModelCriterion.BIC, results_orig.ModelCriterion.LogLikelihood, results_orig.ModelCriterion.Deviance);
            %    results_VTURU.ModelCriterion.AIC, results_VTURU.ModelCriterion.BIC, results_VTURU.ModelCriterion.LogLikelihood, results_VTURU.ModelCriterion.Deviance);


            % RLPFC
            fprintf('\\multicolumn{6}{|c|}{\\textbf{$\\widehat{\\text{RU}}$ and $\\widehat{\\text{TU}}$ from right RLPFC}} \\\\ \\hline\n');

            RU_roi_idx = 1;
            TU_roi_idx = 2;
            load('univariate_decoder_refactored_roiglm-1_badre_glm45_RU_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_both');
            fprintf('Hybrid augmented with $\\widehat{\\text{RU}}$ & V + RU + V/TU + $\\widehat{\\text{RU}}$   &  %.2f  &  %.2f  &  %.2f   &  %.2f     \\\\ \\hline\n', ...
                results_both{RU_roi_idx}.ModelCriterion.AIC, results_both{RU_roi_idx}.ModelCriterion.BIC, results_both{RU_roi_idx}.ModelCriterion.LogLikelihood, results_both{RU_roi_idx}.ModelCriterion.Deviance);

            load('univariate_decoder_refactored_roiglm-1_badre_glm45_TU_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_both');
            fprintf('Hybrid augmented with $\\widehat{\\text{TU}}$ & V + RU + V/TU + V/$\\widehat{\\text{TU}}$   &   %.2f  &  %.2f  &  %.2f   &  %.2f     \\\\ \\hline\n', ...
                results_both{TU_roi_idx}.ModelCriterion.AIC, results_both{TU_roi_idx}.ModelCriterion.BIC, results_both{TU_roi_idx}.ModelCriterion.LogLikelihood, results_both{TU_roi_idx}.ModelCriterion.Deviance);


            % DLPFC
            fprintf('\\multicolumn{6}{|c|}{\\textbf{$\\widehat{\\text{RU}}$ and $\\widehat{\\text{TU}}$ from right DLPFC}} \\\\ \\hline\n');

            load('univariate_decoder_refactored_roiglm-1_dlpfc_glm45_RU_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_both');
            fprintf('Hybrid augmented with $\\widehat{\\text{RU}}$ & V + RU + V/TU + $\\widehat{\\text{RU}}$    &   %.2f  &  %.2f  &  %.2f   &  %.2f       \\\\ \\hline\n', ...
                results_both{RU_roi_idx}.ModelCriterion.AIC, results_both{RU_roi_idx}.ModelCriterion.BIC, results_both{RU_roi_idx}.ModelCriterion.LogLikelihood, results_both{RU_roi_idx}.ModelCriterion.Deviance);

            load('univariate_decoder_refactored_roiglm-1_dlpfc_glm45_TU_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_both');
            fprintf('Hybrid augmented with $\\widehat{\\text{TU}}$ & V + RU + V/TU + V/$\\widehat{\\text{TU}}$  &   %.2f  &  %.2f  &  %.2f   &  %.2f      \\\\ \\hline \\hline\n', ...
                results_both{TU_roi_idx}.ModelCriterion.AIC, results_both{TU_roi_idx}.ModelCriterion.BIC, results_both{TU_roi_idx}.ModelCriterion.LogLikelihood, results_both{TU_roi_idx}.ModelCriterion.Deviance);


            % both
            fprintf('\\multicolumn{6}{|c|}{\\textbf{$\\widehat{\\text{RU}}$ from right RLPFC and $\\widehat{\\text{TU}}$ from right DLPFC}} \\\\ \\hline\n');
            load('univariate_decoder_both_badre_glm45_RUroi=1_TUroi=2_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1.mat', 'results_all');
            fprintf('Hybrid augmented with $\\widehat{\\text{RU}}$ and $\\widehat{\\text{TU}}$ & V + RU + V/TU + $\\widehat{\\text{RU}}$ + V/$\\widehat{\\text{TU}}$   &   %.2f  &  %.2f  &  %.2f   &  %.2f    \\\\ \\hline \\hline\n', ...
                results_all.ModelCriterion.AIC, results_all.ModelCriterion.BIC, results_all.ModelCriterion.LogLikelihood, results_all.ModelCriterion.Deviance);


            % DV
            DV_roi_idx = 1;
            fprintf('\\multicolumn{5}{|c|}{\\textbf{$\\widehat{\\text{DV}}$ from left M1}} \\\\ \\hline\n');
            load('univariate_decoder_refactored_roiglm29_DV_glm29_DV_orth=0_lambda=1.000000_standardize=2_mixed=1_corr=0_extent=100_Num=1_intercept=1_flip=1_doCV=0_gn=0_s=10.0.mat', 'results_both');
            fprintf('Hybrid augmented with $\\widehat{\\text{DV}}$ & V + RU + V/TU + $\\widehat{\\text{DV}}$   &   %.2f  &  %.2f  &  %.2f   &  %.2f    \\\\ \\hline\n', ...
                results_both{DV_roi_idx}.ModelCriterion.AIC, results_both{DV_roi_idx}.ModelCriterion.BIC, results_both{DV_roi_idx}.ModelCriterion.LogLikelihood, results_both{DV_roi_idx}.ModelCriterion.Deviance);




        case 'fig:DV'

            % DV contrast from GLM 29
            %
            figure('pos', [100 100 650 180]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(1,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.60;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(2,1, 1, 'position', pos);

            %PICpng = imread('images/badre_RLPFC.png');   %  <-- to cross-check ; bspmview('masks/badre_rlpfc_36_56_-8_r=10.0mm.nii', '../glmOutput/mean.nii')
            PICpng = imread('images/DV_corr.png'); % extent >= 100

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  


            title('Decision value (corr.)', 'FontSize', fontsize);


            subplot(1,4,3);

            DV_roi_idx = 2;

            load('cross_subj_perf_bic_roiglm29_DV_glm29_sphere_standardize=0_corr=0_extent=100.mat');
            bic = all_bic{DV_roi_idx}';

            % corrs
            %
            [r, p] = corr(perf, bic);
            fprintf('r = %.2f, p = %.3f, Pearson correlation across %d subjects\n', r, p, length(perf));

            scatter(bic, perf);
            hold on;
            lsline;
            text(0, 0.88, sprintf('r = %.2f,\np = %.2f', r, p));
            hold off;

            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XLim',[-8e5 6e5], 'Ylim', [0.64 0.9]);
            xlabel('GLM 2 BIC','FontSize',axisfontsize);
            xticklabels({'-500000', '0', '500000'});
            ylabel('P(better option)','FontSize',axisfontsize);
            %h = get(gca, 'xlabel');
            %get(h,'Position')
            %set(h,'Position',get(h,'Position').*[0 1 1]);



            subplot(1,4,4);

            load('cross_subj_lik_bic_roiglm29_DV_glm29_sphere_standardize=0_corr=0_extent=100.mat');
            bic = all_bic{DV_roi_idx}';

            % corrs
            %
            [r, p] = corr(loglik, bic);
            fprintf('r = %.2f, p = %.3f, Pearson correlation across %d subjects\n', r, p, length(loglik));

            scatter(bic, loglik);
            hold on;
            lsline;
            text(1.5e5, -60, sprintf('r = %.2f,\np = %.2f', r, p));
            hold off;

            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XLim',[-8e5 6e5]);
            xlabel('GLM 2 BIC','FontSize',axisfontsize);
            xticklabels({'-500000', '0', '500000'});
            ylabel('Model log likelihood','FontSize',axisfontsize);


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.08, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.71, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/DV', '-dpdf');








        case 'RU_old' % Figure3
            % RU contrast 
            %
            figure('pos', [100 100 350 800]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(2,1,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.85;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(2,1, 1, 'position', pos);

            %PICpng = imread('images/badre_RLPFC.png');   %  <-- to cross-check 
            %PICpng = imread('images/RU-trial.png');
            %PICpng = imread('images/RU-trial_100.png'); % extent >= 100
            PICpng = imread('images/RU_100.png'); % extent >= 100

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  

            % Badre ROI
            %centx = x * 0.81;
            %centy = y * 0.29;

            % our ROI
            centx = x * 0.79;
            centy = y * 0.30;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            title('RU (uncorr.; extent >= 100)', 'FontSize', fontsize);


            subplot(4,2,5);

            RU_roi_idx = 1;

            %load('main_effect_glm21_RU_RU_-_trial.mat');
            load('main_effect_roiglm36_RU_glm36_RU_corr=0_extent=100_Num=1.mat');
            beta(1) = m(RU_roi_idx);
            ci(1) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(1) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{1} = bs{RU_roi_idx};
            pp(1) = p_uncorr(RU_roi_idx);

            %load('main_effect_glm21_TU_RU_-_trial.mat');
            load('main_effect_roiglm36_RU_glm36_TU_corr=0_extent=100_Num=1.mat');
            beta(2) = m(RU_roi_idx);
            ci(2) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(2) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{2} = bs{RU_roi_idx};
            pp(2) = p_uncorr(RU_roi_idx);
            
            % paired t-test = RU - TU contrast
            [h, p, ci, stats] = ttest(betas{1}, betas{2});
            fprintf('paired t-test RU - TU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            %for s = 1:length(betas{1})
            %    plot([1 2], [betas{1}(s) betas{2}(s)], '-o', 'MarkerSize', 2,'MarkerFaceColor','k', 'Color', 'k');
            %end

            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
            end
            y = max(beta + err + 0.06);
            line([1 2], [y y], 'color', 'black');
            text(1.45, y + 0.02, significance(p));
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.05 0.25]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title({'Main effect', 'RLPFC (R) [34 48 -8]'}, 'FontSize', axisfontsize);

          
            subplot(4,2,6);

            %{
            %load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_roiglm36_RU_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            pp(1) = p_corr(RU_roi_idx);

            %load('univariate_decoder_glm21_TU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_roiglm36_RU_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            pp(2) = p_corr(RU_roi_idx);

            save('Figure3C.mat', 'beta', 'err', 'ci');
            %}
            pp(1) = 0.04; % TODO get from paper; re-run stuff on cluster then recalc properly & save shit... (rename .mat files on cluster)
            pp(2) = 0.9;

            err = [];
            load Figure3C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 5, significance(pp(i)));
            end
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title({'Decoding', 'RLPFC (R) [34 48 -8]'}, 'FontSize', axisfontsize);
            ylim([-20 30]);

            


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/RU_old', '-dpdf');




        case 'TU_old' % Figure4
            % TU contrast
            %

            figure('pos', [100 100 350 800]);
            %figure('pos', [100 100 350 720]);


            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
           
            h = subplot(2,1,1);
            pos = get(h, 'position');
            pos
            %pos(1) = pos(1) * 1.9;
            %pos(2) = pos(2) * 0.92;
            %pos(3) = pos(3) * 1.0 * 2/3;
            %pos(4) = pos(4) * 1.0 * 2/3;
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.85;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;

            subplot(2,1, 1, 'position', pos);
            %PICpng = imread('images/TU-trial.png');
            %PICpng = imread('images/TU-trial_100.png'); % extent >= 100
            PICpng = imread('images/TU_100.png'); % extent >= 100



            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('TU (uncorr.; extent >= 100)', 'FontSize', fontsize);
            hold on;

            xs = [0.09, 0.16];
            ys = [0.18, 0.17];
            colors = {[0.99 0.99 0.99], [0.50 0.99 0.50]};
            for i = 1:length(xs)
                centx = x * xs(i);
                centy = y * ys(i);

                r = 50;
                theta = 0 : (2 * pi / 10000) : (2 * pi);
                pline_x = r * cos(theta) + centx;
                pline_y = r * sin(theta) + centy;
                k = ishold;
                plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', colors{i});
            end
            hold off;



            TU_rois = [2, 8];
            dec_p_RU = [0.2 0.2];
            dec_p_TU = [0.006 0.006];
            ROI_names = {'DLPFC (L) [-42 4 28]', 'DLPFC (L) [-38 36 34]'}
            for ri = 1:length(TU_rois)
                TU_roi_idx = TU_rois(ri); 

                subplot(4,2,5 + 2 * (ri - 1));


                %load('main_effect_glm21_RU_TU_-_trial.mat');
                load('main_effect_roiglm36_TU_glm36_RU_corr=0_extent=100_Num=1.mat');
                beta(1) = m(TU_roi_idx);
                ci(1) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
                err(1) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
                betas{1} = bs{TU_roi_idx};
                pp(1) = p_uncorr(TU_roi_idx);

                %load('main_effect_glm21_TU_TU_-_trial.mat');
                load('main_effect_roiglm36_TU_glm36_TU_corr=0_extent=100_Num=1.mat');
                beta(2) = m(TU_roi_idx);
                ci(2) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
                err(2) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
                betas{2} = bs{TU_roi_idx};
                pp(2) = p_uncorr(TU_roi_idx);

                disp(region{TU_roi_idx});

                % paired t-test = RU - TU contrast
                [h, p, ci, stats] = ttest(betas{2}, betas{1});
                fprintf('paired t-test TU - RU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

                plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
                hold on;
                errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                for i = 1:2
                    text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
                end
                y = max(beta + err + 0.06);
                line([1 2], [y y], 'color', 'black');
                text(1.45, y + 0.02, significance(p));
                hold off;
                set(gca,'TickLabelInterpreter','latex');
                set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
                ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
                title({'Main effect', ROI_names{ri}}, 'FontSize', axisfontsize);
                ylim([-0.2 0.25]);

             
                
                subplot(4,2,6 + 2 * (ri - 1));

                %{
                %load('univariate_decoder_glm21_RU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
                load('univariate_decoder_roiglm36_TU_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
                [b,~,s] = fixedEffects(results_both{TU_roi_idx});
                beta(1) = b(4);
                err(1) = s.SE(4);
                ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

                %load('univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
                load('univariate_decoder_roiglm36_TU_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
                [b,~,s] = fixedEffects(results_both{TU_roi_idx});
                beta(2) = b(4);
                err(2) = s.SE(4);
                ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

                save('Figure4C.mat', 'beta', 'err', 'ci');
                %}
                pp(1) = dec_p_RU(ri); % TODO get from paper; re-run stuff on cluster then recalc properly & save shit... (rename .mat files on cluster)
                pp(2) = dec_p_TU(ri);

                err = [];
                load Figure4C;

                plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
                hold on;
                errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                for i = 1:2
                    text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 2, significance(pp(i)));
                end
                hold off;
                set(gca,'TickLabelInterpreter','latex');
                set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-11 4]);
                ylabel('Regression coefficient (w)','FontSize',axisfontsize);
                title({'Decoding', ROI_names{ri}}, 'FontSize', axisfontsize);
                ylim([-20 5]);

            end

            %{
            subplot(4,3,9);

            load cross_subject_glm21_TU_TU_-_trial_sphere.mat;
            %pos = get(h, 'position');
            %pos(3) = pos(3) * 0.6;
            %pos(1) = pos(1) * 0.9;
            %subplot(3,2, 4, 'position', pos);
            scatter(all_b{6}', w(:,3), 15);
            lsline;
            xlabel('\beta_{TU}', 'FontSize', axisfontsize);
            ylabel('w_3', 'FontSize', axisfontsize);
            t = title({'Variability'}, 'interpreter', 'none', 'FontSize', fontsize);
            set(t,'position',get(t,'position')+[0.03 0 0]);

            str = sprintf('r = %.1f, p = %.2f', r(6), p_uncorr(6));
            text(0.02,0.002, str, 'FontSize', 8);
            %}

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            %text(0.08, 0.78, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            %text(0.08, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            %text(0.39, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');
            %text(0.65, 0.52, 'D', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.78, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.29, 'D', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.29, 'E', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/TU_old', '-dpdf');








        case 'Figure3_old'
            % RU - trial contrast 
            %
            figure('pos', [100 100 650 160]);
            %figure;

            fontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(1,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.93;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(1,2, 1, 'position', pos);

            %PICpng = imread('images/badre_RLPFC.png');   %  <-- to cross-check 
            PICpng = imread('images/RU-trial.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  

            % Badre ROI
            %centx = x * 0.81;
            %centy = y * 0.29;

            % our ROI
            centx = x * 0.80;
            centy = y * 0.30;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            title('RU - trial', 'FontSize', fontsize);


            subplot(1,4,3);

            load('main_effect_glm21_RU_RU_-_trial.mat');
            beta(1) = m(2);
            ci(1) = (cis{2}(2) - cis{2}(1)) / 2;
            err(1) = stat{2}.sd / sqrt(stat{2}.df + 1);

            load('main_effect_glm21_TU_RU_-_trial.mat');
            beta(2) = m(2);
            ci(2) = (cis{2}(2) - cis{2}(1)) / 2;
            err(2) = stat{2}.sd / sqrt(stat{2}.df + 1);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

          
            subplot(1,4,4);

            %{
            load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{2});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            load('univariate_decoder_glm21_TU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{2});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Figure3C.mat', 'beta', 'err', 'ci');
            %}

            err = [];
            load Figure3C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);
            


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Figure3', '-dpdf');




        case 'Figure4_old'
            % TU - trial contrast
            %

            figure('pos', [100 100 650 160]);


            fontsize = 11;
            markersize = 6;
            linewidth = 2;
           
            h = subplot(1,2,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.93;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;

            subplot(1,2, 1, 'position', pos);
            PICpng = imread('images/TU-trial.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('TU - trial', 'FontSize', fontsize);

            centx = x * 0.13;
            centy = y * 0.28;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;


            subplot(1,4,3);

            load('main_effect_glm21_RU_TU_-_trial.mat');
            beta(1) = m(end);
            ci(1) = (cis{end}(2) - cis{end}(1)) / 2;
            err(1) = stat{end}.sd / sqrt(stat{end}.df + 1);

            load('main_effect_glm21_TU_TU_-_trial.mat');
            beta(2) = m(end);
            ci(2) = (cis{end}(2) - cis{end}(1)) / 2;
            err(2) = stat{end}.sd / sqrt(stat{end}.df + 1);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('Insula (L)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

         
            
            subplot(1,4,4);

            %{
            load('univariate_decoder_glm21_RU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{end});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            load('univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{end});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Figure4C.mat', 'beta', 'err', 'ci');
            %}

            err = [];
            load Figure4C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('Insula (L)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-11 4]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Figure4', '-dpdf');


        case 'fig:corr' % FigureS1
            % corrected contrasts
            %

            figure('pos', [100 100 520 450]);


            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
            pos_scale = [1.0 1.0 1.2 1.2];
          
            % RU
            %
            h = subplot(2,2,1);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            pos(2) = pos(2) * 0.8;
            subplot(2,2, 1, 'position', pos);
            PICpng = imread('images/RU_corr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('Relative uncertainty (corr.)', 'FontSize', fontsize);

            % TODO dedupe with Figure3
            % our ROI
            centx = x * 0.79;
            centy = y * 0.30;

            %{
            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;
            %}


            % TU 
            %
            h = subplot(2,2,2);
            pos = get(h, 'position');
            pos(2) = pos(2) * 0.8;
            pos = pos .* pos_scale;
            subplot(2,2, 2, 'position', pos);
            PICpng = imread('images/TU_corr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('Total uncertainty (corr.)', 'FontSize', fontsize);

            % TODO dedupe with Figure4
            hold on;

            %{
            xs = [0.09, 0.16];
            ys = [0.18, 0.17];
            colors = {[0.99 0.99 0.99], [0.50 0.99 0.50]};
            for i = 1:length(xs)
                centx = x * xs(i);
                centy = y * ys(i);

                r = 50;
                theta = 0 : (2 * pi / 10000) : (2 * pi);
                pline_x = r * cos(theta) + centx;
                pline_y = r * sin(theta) + centy;
                k = ishold;
                plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', colors{i});
            end
            hold off;
            %}


            % RU - TU
            %
            %{
            h = subplot(2,2,3);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            subplot(2,2, 3, 'position', pos);
            PICpng = imread('images/RU-TU.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('RU - TU', 'FontSize', fontsize);

            centx = x * 0.80;
            centy = y * 0.30;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;


            % TU - RU
            %
            h = subplot(2,2,4);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            subplot(2,2, 4, 'position', pos);
            PICpng = imread('images/TU-RU.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('TU - RU', 'FontSize', fontsize);

            centx = x * 0.13;
            centy = y * 0.28;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;
            %}

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.85, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.54, 0.85, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            %text(0.13, 0.50, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            %text(0.57, 0.50, 'D', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/corr', '-dpdf');




        case 'badreRU' % FigureS2
            % Badre RLPFC 
            % parallels Figure 3
            %
            
            % EXPT = exploration_expt()
            % struc = fullfile(EXPT.modeldir, 'mean.nii')
            % masks = badre_2012_create_masks(false);
            % bspmview(masks{1}, struc);
            %

            figure('pos', [100 100 350 720]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
          

            h = subplot(2,1,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.85;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;
            subplot(2,1, 1, 'position', pos);

            PICpng = imread('images/badre_RLPFC.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  

            centx = x * 0.81;
            centy = y * 0.29;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            %plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            title('RU (Badre et al. 2012)', 'FontSize', fontsize);

            RU_roi_idx = 1;

            subplot(4,2,5);

            load('main_effect_roiglm-1_badre_glm36_RU_corr=0_extent=100_Num=3.mat');
            beta(1) = m(RU_roi_idx);
            ci(1) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(1) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{1} = bs{RU_roi_idx};
            pp(1) = p_uncorr(RU_roi_idx);

            load('main_effect_roiglm-1_badre_glm36_TU_corr=0_extent=100_Num=3.mat');
            beta(2) = m(RU_roi_idx);
            ci(2) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(2) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{2} = bs{RU_roi_idx};
            pp(2) = p_uncorr(RU_roi_idx);

            % paired t-test = RU - TU contrast
            [h, p, ci, stats] = ttest(betas{1}, betas{2});
            fprintf('paired t-test RU - TU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
            end
            y = max(beta + err + 0.06);
            line([1 2], [y y], 'color', 'black');
            text(1.45, y + 0.02, significance(p));
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.25]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title({'Main effect', 'RLPFC (R) [36 56 -8]'}, 'FontSize', axisfontsize);

          
            subplot(4,2,6);

            %{
            %load('univariate_decoder_glm21_RU_badre_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_RU_badre_norm=4_orth=1_lambda=1.000000.mat');
            load('univariate_decoder_roiglm-1_badre_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=3.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            %load('univariate_decoder_glm21_TU_badre_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_TU_badre_norm=4_orth=1_lambda=1.000000.mat');
            load('univariate_decoder_roiglm-1_badre_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=3.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Supp_Figure1C.mat', 'beta', 'err', 'ci');
            %}
            pp(1) = 0.02; % TODO get from paper; re-run stuff on cluster then recalc properly & save shit... (rename .mat files on cluster)
            pp(2) = 0.9;

            err = [];
            load Supp_Figure1C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 1, significance(pp(i)));
            end
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-7 4]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title({'Decoding', 'RLPFC (R) [36 56 -8]'}, 'FontSize', axisfontsize);
            

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/badreRU', '-dpdf');


            


        case 'badreTU'  % FigureS3
            % Badre DLPFC
            % Parallels Figure 4
            %

            figure('pos', [100 100 350 720]);
            %figure;

            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
           
            h = subplot(2,1,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1;
            pos(2) = pos(2) * 0.85;
            pos(3) = pos(3) * 1.0;
            pos(4) = pos(4) * 1.0;

            subplot(2,1, 1, 'position', pos);
            PICpng = imread('images/badre_DLPFC.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('TU (Badre et al. 2012)', 'FontSize', fontsize);

            centx = x * 0.13;
            centy = y * 0.28;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            %plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;

            TU_roi_idx = 1;

            subplot(4,2,5);

            load('main_effect_roiglm-1_dlpfc_glm36_RU_corr=0_extent=100_Num=3.mat');
            beta(1) = m(TU_roi_idx);
            ci(1) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
            err(1) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
            betas{1} = bs{TU_roi_idx};
            pp(1) = p_uncorr(TU_roi_idx);

            load('main_effect_roiglm-1_dlpfc_glm36_TU_corr=0_extent=100_Num=3.mat');
            beta(2) = m(TU_roi_idx);
            ci(2) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
            err(2) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
            betas{2} = bs{TU_roi_idx};
            pp(2) = p_uncorr(TU_roi_idx);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 0.02, significance(pp(i)));
            end
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.05 0.1]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title({'Main effect', 'DLPFC (R) [40 30 34]'}, 'FontSize', axisfontsize);

         
            
            subplot(4,2,6);

            %{
            %load('univariate_decoder_glm21_RU_dlpfc_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat');
            load('univariate_decoder_roiglm-1_dlpfc_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=3.mat');
            [b,~,s] = fixedEffects(results_both{TU_roi_idx});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            %load('univariate_decoder_glm21_TU_dlpfc_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat');
            load('univariate_decoder_roiglm-1_dlpfc_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=3.mat');
            [b,~,s] = fixedEffects(results_both{TU_roi_idx});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Supp_Figure2C.mat', 'beta', 'err', 'ci');
            %}
            pp(1) = 0.9; % TODO get from paper; re-run stuff on cluster then recalc properly & save shit... (rename .mat files on cluster)
            pp(2) = 0.06;

            err = [];
            load Supp_Figure2C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            for i = 1:2
                text(i - 0.03 * length(significance(pp(i))), beta(i) + err(i) + 2, significance(pp(i)));
            end
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 12]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title({'Decoding', 'DLPFC (R) [40 30 34]'}, 'FontSize', axisfontsize);


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');


            print('images/badreTU', '-dpdf');


        case 'fig:1reg_only'
            % single-regressor GLMs (67 68 69 70) TODO
            %

            figure('pos', [100 100 520 450]);


            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
            pos_scale = [1.0 1.0 1.2 1.2];

            % RU only
            %
            h = subplot(2,2,1);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            pos(2) = pos(2) * 0.8;
            subplot(2,2, 1, 'position', pos);
            PICpng = imread('images/RU_67_uncorr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title({'Relative uncertainty alone', '(uncorr.)'}, 'FontSize', fontsize);

            % TODO dedupe with Figure3
            %{
            centx = x * 0.79;
            centy = y * 0.30;

            r = 50;
            hold on;
            theta = 0 : (2 * pi / 10000) : (2 * pi);
            pline_x = r * cos(theta) + centx;
            pline_y = r * sin(theta) + centy;
            k = ishold;
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            hold off;
            %}


          
            % TU only
            %
            h = subplot(2,2,2);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            pos(2) = pos(2) * 0.8;
            subplot(2,2, 2, 'position', pos);
            PICpng = imread('images/TU_68_uncorr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title({'Total uncertainty alone', '(uncorr.)'}, 'FontSize', fontsize);

            hold on;

            %{
            xs = [0.09, 0.16];
            ys = [0.18, 0.17];
            colors = {[0.99 0.99 0.99], [0.50 0.99 0.50]};
            for i = 1:length(xs)
                centx = x * xs(i);
                centy = y * ys(i);

                r = 50;
                theta = 0 : (2 * pi / 10000) : (2 * pi);
                pline_x = r * cos(theta) + centx;
                pline_y = r * sin(theta) + centy;
                k = ishold;
                plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', colors{i});
            end
            hold off;
            %}



            % V only
            %
            h = subplot(2,2,3);
            pos = get(h, 'position');
            pos(2) = pos(2) * 0.8;
            pos = pos .* pos_scale;
            subplot(2,2, 3, 'position', pos);
            PICpng = imread('images/V_69_uncorr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title({'Value difference alone', ' (uncorr.)'}, 'FontSize', fontsize);


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.88, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.55, 0.88, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.10, 0.50, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            %text(0.57, 0.50, 'D', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/1reg_only', '-dpdf');





        case 'Supp_Figure3_old'
            % Cross-subject
            %


            figure('pos', [100 100 300 250]);


            fontsize = 12;

            %{
            h = subplot(1,2,1);
            load cross_subj_vbm_glm21_RU_RU_-_trial_sphere_standardize=0.mat;
            %pos = get(h, 'position');
            %pos(3) = pos(3) * 0.6;
            %pos(1) = pos(1) * 0.9;
            %subplot(3,2, 4, 'position', pos);
            scatter(all_b{2}', w(:,2));
            lsline;
            xlabel('Grey matter density');
            ylabel('w_2');
            title('RLPFC (R)', 'interpreter', 'none', 'FontSize', fontsize);

            str = sprintf('r = %.2f, p = %.3f', r(2), p_uncorr(2));
            text(1.15,0.004, str, 'FontSize', 9);

            %}




            %h = subplot(1,2,2);

            load cross_subject_glm21_TU_TU_-_trial_sphere.mat;
            %pos = get(h, 'position');
            %pos(3) = pos(3) * 0.6;
            %pos(1) = pos(1) * 0.9;
            %subplot(3,2, 4, 'position', pos);
            scatter(all_b{6}', w(:,3));
            lsline;
            xlabel('\beta_{TU}');
            ylabel('w_3');
            title('Insula (L)', 'interpreter', 'none', 'FontSize', fontsize);

            str = sprintf('r = %.2f, p = %.3f', r(6), p_uncorr(6));
            text(-0.13,0.004, str, 'FontSize', 9);

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            %text(0.13, 0.65, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            %text(0.48, 0.65, 'B', 'FontSize', 20, 'FontWeight', 'bold');


            print('images/Supp_Figure3', '-dpdf');




        case 'TableS2'
            % RU 

            %ccnl_view(exploration_expt(), 21, 'RU - trial');
            ccnl_results_table('AAL2', 'peak', exploration_expt(), 36, 'RU', 0.001, '+/-', 0.05, 20, 1, false, 100); % extent >= 100

        case 'TableS3'
            % TU trial

            tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 36, 'TU', 0.001, '+/-', 0.05, 20, 1, false, 100); % extent >= 100

        case 'TableS4'
            % RU univaraite decoder

            %{
            %load univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat;
            load('univariate_decoder_roiglm36_RU_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            save('Supp_Table1.mat', 'p_uncorr', 'p_corr', 'comp');
            %}

            load('Supp_Table1.mat');

            ROI_idx = 1;

            fprintf('\n\n\n');
            fprintf('RLPFC (R) p-value = %.4f, loglik = -3748.3 (glm) vs. -3745.8 (altglm); BICs = 7524 vs. 7528.1; AICs = 7502.5 vs. 7499.5\n', p_uncorr(ROI_idx));
            fprintf('\n\n\n');

            fprintf('Hybrid model       &  %d  &  %.2f  &  %.2f  &  %.2f  &         &  \\\\ \n', comp{ROI_idx}.DF(1), comp{ROI_idx}.AIC(1), comp{ROI_idx}.BIC(1), comp{ROI_idx}.LogLik(1));
            fprintf('Augmented RU model &  %d  &  %.2f  &  %.2f  &  %.2f  &   %.2f  &  p = %.3f     \\\\ \n', comp{ROI_idx}.DF(2), comp{ROI_idx}.AIC(2), comp{ROI_idx}.BIC(2), comp{ROI_idx}.LogLik(2), comp{ROI_idx}.LRStat(2), comp{ROI_idx}.pValue(2));



        case 'TableS5'
            % TU with univariate decoder p's

            %ccnl_view(exploration_expt(), 21, 'TU - trial');
            %tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+', 0.05, 20, 1);
            tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 36, 'TU', 0.001, '+', 0.05, 20, 1, false, 100);

            %{
            %load univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat;
            load('univariate_decoder_roiglm36_TU_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            save('Table2.mat', 'p_uncorr', 'p_corr', 'comp');
            %}

            % uncorrected
            % load univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0.mat

            %load('Table2.mat');
            %load('Table2_uncorr.mat'); % uncorrected
            load('Supp_Table2.mat');


            table(p_uncorr, p_corr)
            fprintf('\n\n\n');

            %extent = 100; % manually apply correction with custom extent
            %which = [tab{:,4}] >= extent;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);

            % \textbf{Brain region} & \textbf{MNI coord.} & \textbf{AIC} & \textbf{BIC}  & \textbf{LL}  & \textbf{LR-stat} & \textbf{p-value (uncorr.)} & \textbf{p-value (corr.)} \\ \hline
            for i = 1:length(p_uncorr)
                %if ~which(i)
                %    continue;
                %end
                fprintf('%s & %d %d %d & %.2f & %.2f & %.2f & %.2f &', tab{i,2}, tab{i,end-2}, tab{i,end-1}, tab{i,end}, comp{i}.AIC(2), comp{i}.BIC(2), comp{i}.LogLik(2), comp{i}.LRStat(2));

                p = p_uncorr(i);
                if p > 0.0001
                    p_string1 = sprintf('p = %.4f', p);
                else
                    p_string1 = sprintf('p < 10^{%.0f}', ceil(log10(p)));
                end


                p = p_corr(i);
                if p > 0.0001
                    p_string = sprintf('p = %.4f', p);
                else
                    p_string = sprintf('p < 10^{%.0f}', ceil(log10(p)));
                end

                if p < 0.05
                    p_string = ['*', p_string];
                    p_string1 = ['*', p_string1];
                end
                fprintf('$%s$ & $%s$ \\\\ \\hline \n', p_string1, p_string);
            end


        case 'Table6'
            % Cross-subject analysis for TU

            %tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+/-', 0.05, 20, 1);
            tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+', 0.05, 20, 1, false); % uncorrected

            %load('Table2.mat');
            load('Table2_uncorr.mat'); % uncorrected

            % manually apply correction with custom extent
            extent = 100;
            which = [tab{:,4}] >= extent;
            p_corr = 1 - (1 - p_uncorr) .^ sum(which);

            idx = which' & (p_corr < 0.05); % which were significant for univariate decoder

            load cross_subject_glm21_TU_TU_-_trial_sphere.mat;
            load('cross_subject_glm21_TU_TU_-_trial_sphere_standardize=2_corr=0.mat'); % uncorr

            save shit.mat

            names = region(idx);
            names = tab(idx,2);
            p_uncorr = ps(idx);
            p_corr = 1 - (1 - p_uncorr) .^ length(p_uncorr);
            r = rs(idx);
            mni = mni(idx,:);

            table(names, mni, r, p_uncorr, p_corr)

            for i = 1:length(names)
                if p_corr(i) < 0.05
                    star = '*';
                else
                    star = '';
                end
                fprintf('%s & %d %d %d & %.2f & %sp = %.3f & %sp = %.3f \\\\ \\hline \n', ...
                    names{i}, ...
                    mni(i,1), mni(i,2), mni(i,3), ...
                    r(i), ...
                    star, p_uncorr(i), ...
                    star, p_corr(i));
                end



        case 'TableS7'

            % final augmented model TODO finish
            load univariate_decoder_both_glm21_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100.mat;





        case 'Sam_Figure1'
            fontsize = 25;
            linewidth = 5;

            x = linspace(-10,10,100);
            plot(x,normpdf(x,2,4),'-k','LineWidth',linewidth);
            hold on;
            x = -1; y = normpdf(x,x,4);
            plot([x x],[0 y],'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5])
            text(2-0.4,y+0.01,'R','FontSize',fontsize);
            text(-1-0.4,y+0.01,'S','FontSize',fontsize);
            set(gca,'FontSize',fontsize,'YLim',[0 0.12]);
            ylabel('Probability density','FontSize',fontsize);
            xlabel('Reward','FontSize',fontsize);


        case 'Sam_Figure2'
            fontsize = 25;
            linewidth = 4;
            markersize = 8;

            v = linspace(-30,30,8);
            x = v(1:end-1) + diff(v)/2;
           
            data = load_data;
            tbl = data2table(data,1);
            load results_glme_fig3.mat;
            CC = predict(results_VTURU,tbl);
            
            for s = 1:length(data)
                latents = kalman_filter(data(s));
                V = latents.m(~data(s).timeout,1) - latents.m(~data(s).timeout,2);
                C = double(data(s).choice(~data(s).timeout)==1);
                c = CC(tbl.S==s);
                
                for j = 1:length(v)-1
                    for cond = 1:4
                        ix = V>v(j) & V<v(j+1) & data(s).cond(~data(s).timeout)==cond;
                        if ~any(ix)
                            pc0(s,j,cond) = nan;
                            pc1(s,j,cond) = nan;
                        else
                            pc0(s,j,cond) = nanmean(C(ix));
                            pc1(s,j,cond) = nanmean(c(ix));
                        end
                    end
                end
            end
            
            mu = linspace(-3,3,100);
            d = 1;
            p(1,:) = 1-normcdf(0,mu+d,1);
            p(2,:) = 1-normcdf(0,mu,1+d);
            
            figure;
            T = {'Relative uncertainty: intercept shift' 'Total uncertainty: slope shift'};
            for i = 1:2
                subplot(2,2,i);
                plot(mu,p(i,:),'-k','LineWidth',linewidth); hold on;
                plot(mu,1-normcdf(0,mu,1),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                if i==1
                    legend({'RS' 'SR'},'FontSize',fontsize,'Location','East');
                else
                    legend({'RR' 'SS'},'FontSize',fontsize,'Location','East');
                end
                set(gca,'FontSize',fontsize,'XLim',[min(mu) max(mu)],'YLim',[-0.05 1.05]);
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference (V)','FontSize',fontsize);
                title(T{i},'FontSize',fontsize','FontWeight','Bold');
            end
            
            [se,mu] = wse(pc0);
            mu_c = squeeze(nanmean(pc1));
            
            for i = 1:2
                subplot(2,2,i+2);
                
                if i==1
                    p1 = errorbar(x/10,mu(:,1),se(:,1),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,2),se(:,2),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,1),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,2),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RS' 'SR'},'FontSize',fontsize,'Location','East');
                else
                    p1 = errorbar(x/10,mu(:,3),se(:,3),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,4),se(:,4),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,3),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,4),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RR' 'SS'},'FontSize',fontsize,'Location','East');
                end
                
                set(gca,'FontSize',fontsize,'XLim',[-3 3],'YLim',[-0.05 1.05]);
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference','FontSize',fontsize);
            end
            
            
            set(gcf,'Position',[200 200 1000 700]);
            
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.67, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.35, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            

        case 'Sam_Figure4'

            % Probit analysis of conditions
            data = load_data;
            markersize = 12;
            fontsize = 25;
            
            % fit generalized linear mixed effects model
            if ~exist('results_glme_fig4.mat', 'file')
                tbl = data2table(data,1);
                formula = 'C ~ -1 + cond + cond:V + (-1 + cond + cond:V|S)';
                % note: Laplace method not used here because it doesn't seem to complete
                results = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','DummyVarCoding','Full');
                save results_glme_fig4 results
            else
                load results_glme_fig4.mat
            end
            
            % hypothesis tests
            H = [1 -1 0 0 0 0 0 0; 0 0 1 -1 0 0 0 0; 0 0 0 0 1 -1 0 0; 0 0 0 0 0 0 1 -1];   % contrast matrix
            for i=1:4; [p(i),F(i),DF1(i),DF2(i)] = coefTest(results,H(i,:)); end
            disp(['intercept, RS vs. SR: F(',num2str(DF1(1)),',',num2str(DF2(1)),') = ',num2str(F(1)),', p = ',num2str(p(1))]);
            disp(['intercept, RR vs. SS: F(',num2str(DF1(2)),',',num2str(DF2(2)),') = ',num2str(F(2)),', p = ',num2str(p(2))]);
            disp(['slope, RS vs. SR: F(',num2str(DF1(3)),',',num2str(DF2(3)),') = ',num2str(F(3)),', p = ',num2str(p(3))]);
            disp(['slope, RR vs. SS: F(',num2str(DF1(4)),',',num2str(DF2(4)),') = ',num2str(F(4)),', p = ',num2str(p(4))]);

            % plot results
            figure;
            [beta,~,stats] = fixedEffects(results);
            subplot(1,2,1);
            errorbar(beta(1:4),stats.SE(1:4),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5],'YLim', [-1 1]);
            ylabel('Intercept','FontSize',fontsize);
            subplot(1,2,2);
            errorbar(beta(5:8),stats.SE(5:8),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [0 3]);
            ylabel('Slope','FontSize',fontsize);
            set(gcf,'Position',[200 200 1000 400])

        case 'Sam_Figure3'
            
            % Probit analysis of computational variables
            
            load results_glme_fig3
            results = results_VTURU;
                        
            % plot results
            [beta,~,stats] = fixedEffects(results);
            errorbar(beta([3 1 2]),stats.SE([3 1 2]),'ok','MarkerSize',12,'MarkerFaceColor','k');
            set(gca,'FontSize',25,'XTickLabel',{'V' 'RU' 'V/TU'},'XLim',[0.5 3.5],'YLim',[0 3]);
            ylabel('Regression coefficient','FontSize',25);
            set(gcf,'Position',[200 200 500 400])

        case 'Sam_Figure1'
            fontsize = 25;
            linewidth = 5;

            x = linspace(-10,10,100);
            plot(x,normpdf(x,2,4),'-k','LineWidth',linewidth);
            hold on;
            x = -1; y = normpdf(x,x,4);
            plot([x x],[0 y],'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5])
            text(2-0.4,y+0.01,'R','FontSize',fontsize);
            text(-1-0.4,y+0.01,'S','FontSize',fontsize);
            set(gca,'FontSize',fontsize,'YLim',[0 0.12]);
            ylabel('Probability density','FontSize',fontsize);
            xlabel('Reward','FontSize',fontsize);


        case 'Sam_Figure2'
            fontsize = 25;
            linewidth = 4;
            markersize = 8;

            v = linspace(-30,30,8);
            x = v(1:end-1) + diff(v)/2;
           
            data = load_data;
            tbl = data2table(data,1);
            load results_glme_fig3.mat;
            CC = predict(results_VTURU,tbl);
            
            for s = 1:length(data)
                latents = kalman_filter(data(s));
                V = latents.m(~data(s).timeout,1) - latents.m(~data(s).timeout,2);
                C = double(data(s).choice(~data(s).timeout)==1);
                c = CC(tbl.S==s);
                
                for j = 1:length(v)-1
                    for cond = 1:4
                        ix = V>v(j) & V<v(j+1) & data(s).cond(~data(s).timeout)==cond;
                        if ~any(ix)
                            pc0(s,j,cond) = nan;
                            pc1(s,j,cond) = nan;
                        else
                            pc0(s,j,cond) = nanmean(C(ix));
                            pc1(s,j,cond) = nanmean(c(ix));
                        end
                    end
                end
            end
            
            mu = linspace(-3,3,100);
            d = 1;
            p(1,:) = 1-normcdf(0,mu+d,1);
            p(2,:) = 1-normcdf(0,mu,1+d);
            
            figure;
            T = {'Relative uncertainty: intercept shift' 'Total uncertainty: slope shift'};
            for i = 1:2
                subplot(2,2,i);
                plot(mu,p(i,:),'-k','LineWidth',linewidth); hold on;
                plot(mu,1-normcdf(0,mu,1),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                if i==1
                    legend({'RS' 'SR'},'FontSize',fontsize,'Location','East');
                else
                    legend({'RR' 'SS'},'FontSize',fontsize,'Location','East');
                end
                set(gca,'FontSize',fontsize,'XLim',[min(mu) max(mu)],'YLim',[-0.05 1.05]);
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference (V)','FontSize',fontsize);
                title(T{i},'FontSize',fontsize','FontWeight','Bold');
            end
            
            [se,mu] = wse(pc0);
            mu_c = squeeze(nanmean(pc1));
            
            for i = 1:2
                subplot(2,2,i+2);
                
                if i==1
                    p1 = errorbar(x/10,mu(:,1),se(:,1),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,2),se(:,2),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,1),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,2),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RS' 'SR'},'FontSize',fontsize,'Location','East');
                else
                    p1 = errorbar(x/10,mu(:,3),se(:,3),'ok','LineWidth',linewidth,'MarkerFaceColor','k','MarkerSize',markersize); hold on;
                    p2 = errorbar(x/10,mu(:,4),se(:,4),'o','LineWidth',linewidth,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',markersize);
                    plot(x/10,mu_c(:,3),'-k','LineWidth',linewidth); hold on;
                    plot(x/10,mu_c(:,4),'-','LineWidth',linewidth,'Color',[0.5 0.5 0.5]);
                    %legend([p1 p2],{'RR' 'SS'},'FontSize',fontsize,'Location','East');
                end
                
                set(gca,'FontSize',fontsize,'XLim',[-3 3],'YLim',[-0.05 1.05]);
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference','FontSize',fontsize);
            end
            
            
            set(gcf,'Position',[200 200 1000 700]);

    end

end




function draw_circle(centx, centy, r, color)
    theta = 0 : (2 * pi / 10000) : (2 * pi);
    pline_x = r * cos(theta) + centx;
    pline_y = r * sin(theta) + centy;
    k = ishold;
    plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', color);
end
