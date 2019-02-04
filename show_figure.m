function show_figure(fig)

    switch fig
            
        case 'badre_ROI'
            % show the Badre ROI 
            %
            EXPT = exploration_expt();
            struc = fullfile(EXPT.modeldir,'mean.nii');
            masks = badre_2012_create_masks(false);
            bspmview(masks{1}, struc);


        case 'vifs'

            %{
            EXPT = exploration_expt();
            [vifs, names] = ccnl_vifs(EXPT, 36);

            save vifs.mat;
            %}

            load vifs.mat;

            figure;

            % RU
            %
            for s = 1:length(vifs)
                RU_vifs = vifs{s}(contains(names{s}, 'xRU'));

                subplot(8, 4, s);
                plot(RU_vifs, 'ko', 'MarkerFaceColor', [1 0.5 0]);
                hold on;
                plot([0 length(RU_vifs)], [1 1], 'k');
                %plot([0 length(RU_vifs)], [2 2], 'b--');
                %plot([0 length(RU_vifs)], [4 4], 'r--');
                %plot([0 length(RU_vifs)], [8 8], 'r-');
                plot([0 length(RU_vifs)], [10 10], 'r-');
                hold off;

                if s == 2
                    title('RU VIFs for each subject');
                end
                if s == 30
                    xlabel('run');
                end
            end


            % TU
            %
            for s = 1:length(vifs)
                TU_vifs = vifs{s}(contains(names{s}, 'xTU'));

                subplot(8, 4, s);
                plot(TU_vifs, 'ko', 'MarkerFaceColor', [1 0.5 0]);
                hold on;
                plot([0 length(TU_vifs)], [1 1], 'k');
                %plot([0 length(TU_vifs)], [2 2], 'b--');
                %plot([0 length(TU_vifs)], [4 4], 'r--');
                %plot([0 length(TU_vifs)], [8 8], 'r-');
                plot([0 length(TU_vifs)], [10 10], 'r-');
                hold off;

                if s == 2
                    title('TU VIFs for each subject');
                end
                if s == 30
                    xlabel('run');
                end
            end


        case 'recovery'

            figure('pos', [10 10 520 450]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;


            %load recovery.mat % expo [0, infty]
            load recovery_mvnrnd.mat % mvnrnd

            for i = 1:3
                [r,p] = corr(w_orig(:,i), w_rec(:,i));
                fprintf('recovery w_%d: r = %.4f, p = %e\n', i, r, p);

                subplot(2,3,i);
                scatter(w_orig(:,i), w_rec(:,i));
                lsline;
                title(sprintf('w_%d', i));
                xlabel('simulated');
                ylabel('fitted');
            end

            k = 0;
            for i = 1:3
                for j = i+1:3
                    k = k + 1;

                    [r,p] = corr(w_rec(:,i), w_rec(:,j));
                    fprintf('recovery w_%d vs w_%d: r = %.4f, p = %f\n', i, j, r, p);

                    subplot(2,3,k + 3);
                    scatter(w_rec(:,i), w_rec(:,j));
                    lsline;
                    title(sprintf('w_%d vs. w_%d', i, j));
                    xlabel(sprintf('fitted w_%d', i));
                    ylabel(sprintf('fitted w_%d', j));
                end
            end



        case 'learning'

            figure('pos', [10 10 520 450]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;


            conds = {'RS', 'SR', 'RR', 'SS'};

            captions = {'human', 'model'};

            %  learning curves
            %
            data = load_data;
            tbl = data2table(data,1,1); % do standardize! that's how we analyze behavior
            load results_glme_fig3.mat;
            results = results_VTURU;
            y = predict(results, tbl);

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
                            C = double((y(tbl.S == s) > 0.5) == better); % model choices
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
                legend(conds,'FontSize',fontsize,'Location','East');
                %set(gca,'FontSize',fontsize,'XLim',[min(v) max(v)],'YLim',[0 1]);
                ylabel('P(better option)','FontSize',fontsize);
                xlabel('trial','FontSize',fontsize);
                title(captions{human_or_model});
            end


            tbl = data2table(data,1,1);
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



        case 'psycho'

            figure('pos', [10 10 520 450]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;


            conds = {[1 2], [3 4]};
            legends = {{'RS', 'SR'}, {'RR', 'SS'}};

            captions = {'human', 'model'};

            %  psychometric curves
            %

            for human_or_model = 1:2
                for curve = 1:2

                    clear pc;
                
                    data = load_data;
                    tbl = data2table(data,1,1);
                    load results_glme_fig3.mat;
                    results = results_VTURU;
                    y = predict(results, tbl);

                    latents = kalman_filter(data(1));
                    v = linspace(min(latents(1).m(:)),max(latents(1).m(:)),11)';
                    %v = linspace(-25, 25, 11)';

                    %b = [];
                    for s = 1:length(data)
                        latents = kalman_filter(data(s));

                        which = ~data(s).timeout;
                        V = latents.m(which,1) - latents.m(which,2);

                        if human_or_model == 1
                            C = double(data(s).choice(which)==1); % human choices
                        else
                            C = double(y(tbl.S == s) > 0.5); % model choices
                        end

                        
                        for j = 1:length(v)-1
                            ix = V>v(j) & V<v(j+1) & data(s).cond(which) == conds{curve}(1);
                            if ~any(ix)
                                pc(s,j,1) = nan;
                            else
                                pc(s,j,1) = nanmean(C(ix));
                            end
                            
                            ix = V>v(j) & V<v(j+1) & data(s).cond(which) == conds{curve}(2);
                            if ~any(ix)
                                pc(s,j,2) = nan;
                            else
                                pc(s,j,2) = nanmean(C(ix));
                            end
                        end
                    end

                    subplot(2,2, 1 + human_or_model - 1 + 2 * (curve - 1));

                    [se,mu] = wse(pc);
                    x = v(1:end-1) + diff(v)/2;
                    errorbar(x,mu(:,1),se(:,1),'-ok','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor','k'); hold on
                    errorbar(x,mu(:,2),se(:,2),'-o','LineWidth',linewidth,'MarkerSize',markersize,'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);
                    legend(legends{curve},'FontSize',fontsize,'Location','East');
                    set(gca,'FontSize',fontsize,'XLim',[min(v) max(v)],'YLim',[0 1]);
                    ylabel('Choice probability','FontSize',fontsize);
                    xlabel('Expected value difference','FontSize',fontsize);
                    if curve == 1
                        title(captions{human_or_model});
                    end
                end
            end


        case 'Figure1'
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
                ylabel('P(choose left)','FontSize',fontsize);
                xlabel({'Expected value difference,', 'V = Q(left) - Q(right)'},'FontSize',fontsize);
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

            print('images/Figure1', '-dpdf');


        case 'Figure2'
            figure('pos', [10 10 700 200]);

            fontsize = 12;
            linewidth = 3;
            markersize = 6;

            %% -------------------- Sam Figure 4

            % Probit analysis of conditions
            data = load_data;
            
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
            ylim([1 2.5]);


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
            
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.07, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.35, 0.95, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.64, 0.95, 'C', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/Figure2', '-dpdf');









        case 'Figure3'
            % RU contrast 
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

            %PICpng = imread('images/badre_RLPFC.png');   %  <-- to cross-check 
            %PICpng = imread('images/RU-trial.png');
            %PICpng = imread('images/RU-trial_100.png'); % extent >= 100
            PICpng = imread('images/RU.png'); % extent >= 100

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

            title('RU (uncorr.)', 'FontSize', fontsize);


            subplot(4,2,5);

            RU_roi_idx = 1;

            %load('main_effect_glm21_RU_RU_-_trial.mat');
            load('main_effect_roiglm36_RU_glm36_RU_corr=0_extent=100_Num=1.mat');
            beta(1) = m(RU_roi_idx);
            ci(1) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(1) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{1} = bs{RU_roi_idx};

            %load('main_effect_glm21_TU_RU_-_trial.mat');
            load('main_effect_roiglm36_RU_glm36_TU_corr=0_extent=100_Num=1.mat');
            beta(2) = m(RU_roi_idx);
            ci(2) = (cis{RU_roi_idx}(2) - cis{RU_roi_idx}(1)) / 2;
            err(2) = stat{RU_roi_idx}.sd / sqrt(stat{RU_roi_idx}.df + 1);
            betas{2} = bs{RU_roi_idx};
            
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
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title('Main effect', 'FontSize', fontsize);

          
            subplot(4,2,6);

            %{
            %load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_roiglm36_RU_glm36_RU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            %load('univariate_decoder_glm21_TU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_roiglm36_RU_glm36_TU_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0_extent=100_Num=1.mat');
            [b,~,s] = fixedEffects(results_both{RU_roi_idx});
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
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title('Decoding', 'FontSize', fontsize);
            ylim([-20 25]);
            


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Figure3', '-dpdf');




        case 'Figure4'
            % TU contrast
            %

            figure('pos', [100 100 550 800]);


            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
           
            h = subplot(2,1,1);
            pos = get(h, 'position');
            pos
            pos(1) = pos(1) * 1.9;
            pos(2) = pos(2) * 0.92;
            pos(3) = pos(3) * 1.0 * 2/3;
            pos(4) = pos(4) * 1.0 * 2/3;

            subplot(2,1, 1, 'position', pos);
            %PICpng = imread('images/TU-trial.png');
            %PICpng = imread('images/TU-trial_100.png'); % extent >= 100
            PICpng = imread('images/TU.png'); % extent >= 100



            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('TU (uncorr.)', 'FontSize', fontsize);
            hold on;

            xs = [0.09, 0.16];
            ys = [0.18, 0.17];
            for i = 1:length(xs)
                centx = x * xs(i);
                centy = y * ys(i);

                r = 50;
                theta = 0 : (2 * pi / 10000) : (2 * pi);
                pline_x = r * cos(theta) + centx;
                pline_y = r * sin(theta) + centy;
                k = ishold;
                plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [0.99 0.99 0.99]);
            end
            hold off;



            TU_rois = [2, 8];
            for ri = 1:length(TU_rois)
                TU_roi_idx = TU_rois(ri); 

                subplot(4,2,5 + 2 * (ri - 1));

                %TU_roi_idx = 2; % TODO determine & re-render

                %load('main_effect_glm21_RU_TU_-_trial.mat');
                load('main_effect_roiglm36_TU_glm36_RU_corr=0_extent=100_Num=1.mat');
                beta(1) = m(TU_roi_idx);
                ci(1) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
                err(1) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
                betas{1} = bs{TU_roi_idx};

                %load('main_effect_glm21_TU_TU_-_trial.mat');
                load('main_effect_roiglm36_TU_glm36_TU_corr=0_extent=100_Num=1.mat');
                beta(2) = m(TU_roi_idx);
                ci(2) = (cis{TU_roi_idx}(2) - cis{TU_roi_idx}(1)) / 2;
                err(2) = stat{TU_roi_idx}.sd / sqrt(stat{TU_roi_idx}.df + 1);
                betas{2} = bs{TU_roi_idx};

                disp(region{TU_roi_idx});

                plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
                hold on;
                errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                hold off;
                set(gca,'TickLabelInterpreter','latex');
                set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
                ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
                title('Main effect', 'FontSize', fontsize);
                ylim([-0.2 0.2]);

                % paired t-test = RU - TU contrast
                [h, p, ci, stats] = ttest(betas{2}, betas{1});
                fprintf('paired t-test TU - RU: t(%d) = %.4f, p = %.6f\n', stats.df, stats.tstat, p);
             
                
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

                err = [];
                load Figure4C;

                plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
                hold on;
                errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
                hold off;
                set(gca,'TickLabelInterpreter','latex');
                set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-11 4]);
                ylabel('Regression coefficient (w)','FontSize',axisfontsize);
                title('Decoding', 'FontSize', fontsize);
                ylim([-20 5]);

            end
            % TODO prettify

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

            print('images/Figure4', '-dpdf');








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


        case 'FigureS1'
            % corrected contrasts
            %

            figure('pos', [100 100 520 450]);


            fontsize = 14;
            axisfontsize = 11;
            markersize = 6;
            linewidth = 2;
            pos_scale = [1.0 1.0 1.2 1.2];
          
            % RU - trial
            %
            h = subplot(2,2,1);
            pos = get(h, 'position');
            pos = pos .* pos_scale;
            pos(2) = pos(2) * 0.8;
            subplot(2,2, 1, 'position', pos);
            PICpng = imread('images/RU-trial_corr.png');

            [rows columns numberOfColorChannels] = size(PICpng);
            x = columns;
            y = rows;
            imshow(PICpng, 'InitialMagnification', 'fit');  
            title('RU - trial', 'FontSize', fontsize);

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


            % TU - trial
            %
            h = subplot(2,2,2);
            pos = get(h, 'position');
            pos(2) = pos(2) * 0.8;
            pos = pos .* pos_scale;
            subplot(2,2, 2, 'position', pos);
            PICpng = imread('images/TU-trial_corr.png');

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


            % RU - TU
            %
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

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.13, 0.85, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.57, 0.85, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.13, 0.50, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.57, 0.50, 'D', 'FontSize', 25, 'FontWeight', 'bold');

            print('images/FigureS1', '-dpdf');




        case 'FigureS2'
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


            subplot(4,2,5);

            load('main_effect_glm21_RU_badre.mat');
            beta(1) = m(1);
            ci(1) = (cis{1}(2) - cis{1}(1)) / 2;
            err(1) = stat{1}.sd / sqrt(stat{1}.df + 1);

            load('main_effect_glm21_TU_badre.mat');
            beta(2) = m(1);
            ci(2) = (cis{1}(2) - cis{1}(1)) / 2;
            err(2) = stat{1}.sd / sqrt(stat{1}.df + 1);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title('Main effect', 'FontSize', fontsize);

          
            subplot(4,2,6);

            %{
            %load('univariate_decoder_glm21_RU_badre_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_glm21_RU_badre_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{1});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            %load('univariate_decoder_glm21_TU_badre_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            load('univariate_decoder_glm21_TU_badre_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{1});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Supp_Figure1C.mat', 'beta', 'err', 'ci');
            %}

            err = [];
            load Supp_Figure1C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 16]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title('Decoding', 'FontSize', fontsize);
            

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/FigureS2', '-dpdf');


            


        case 'FigureS3'
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


            subplot(4,2,5);

            load('main_effect_glm21_RU_dlpfc.mat');
            beta(1) = m(1);
            ci(1) = (cis{1}(2) - cis{1}(1)) / 2;
            err(1) = stat{1}.sd / sqrt(stat{1}.df + 1);

            load('main_effect_glm21_TU_dlpfc.mat');
            beta(2) = m(1);
            ci(2) = (cis{1}(2) - cis{1}(1)) / 2;
            err(2) = stat{1}.sd / sqrt(stat{1}.df + 1);

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',axisfontsize);
            title('Main effect', 'FontSize', fontsize);

         
            
            subplot(4,2,6);

            %{
            load('univariate_decoder_glm21_RU_dlpfc_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{1});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            load('univariate_decoder_glm21_TU_dlpfc_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat');
            [b,~,s] = fixedEffects(results_both{1});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Supp_Figure2C.mat', 'beta', 'err', 'ci');
            %}

            err = [];
            load Supp_Figure2C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            %errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',axisfontsize,'XTick', [1 2],'XTickLabel',{'$\widehat{RU}$', '$V/\widehat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',axisfontsize);
            title('Decoding', 'FontSize', fontsize);


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.04, 0.80, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.04, 0.52, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.49, 0.52, 'C', 'FontSize', 20, 'FontWeight', 'bold');


            print('images/FigureS3', '-dpdf');





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
            % RU - trial

            %ccnl_view(exploration_expt(), 21, 'RU - trial');
            ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'RU - trial', 0.001, '+/-', 0.05, 20, 1, false, 100); % extent >= 100

        case 'TableS3'
            % TU - trial

            tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+/-', 0.05, 20, 1, false, 100); % extent >= 100

        case 'TableS4'
            % RU univaraite decoder

            %{
            load univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat;
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            save('Supp_Table1.mat', 'p_uncorr', 'p_corr', 'comp');
            %}

            load('Supp_Table1.mat');

            fprintf('\n\n\n');
            fprintf('RLPFC (R) p-value = %.4f, loglik = -3748.3 (glm) vs. -3745.8 (altglm); BICs = 7524 vs. 7528.1; AICs = 7502.5 vs. 7499.5\n', p_uncorr(2));
            fprintf('\n\n\n');

            fprintf('Hybrid model       &  %d  &  %.2f  &  %.2f  &  %.2f  &         &  \\\\ \n', comp{2}.DF(1), comp{2}.AIC(1), comp{2}.BIC(1), comp{2}.LogLik(1));
            fprintf('Augmented RU model &  %d  &  %.2f  &  %.2f  &  %.2f  &   %.2f  &  p = %.3f     \\\\ \n', comp{2}.DF(2), comp{2}.AIC(2), comp{2}.BIC(2), comp{2}.LogLik(2), comp{2}.LRStat(2), comp{2}.pValue(2));



        case 'TableS5'
            % TU - trial with univariate decoder p's

            %ccnl_view(exploration_expt(), 21, 'TU - trial');
            %tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+', 0.05, 20, 1);
            tab = ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+', 0.05, 20, 1, false); % uncorrected

            %{
            load univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=0_mixed=0.mat;
            p_uncorr = p_comp;
            p_corr = 1 - (1 - p_uncorr) .^ numel(p_uncorr);
            save('Table2.mat', 'p_uncorr', 'p_corr', 'comp');
            %}

            % uncorrected
            % load univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0_corr=0.mat

            %load('Table2.mat');
            load('Table2_uncorr.mat'); % uncorrected

            table(p_uncorr, p_corr)
            fprintf('\n\n\n');

            extent = 100; % manually apply correction with custom extent
            which = [tab{:,4}] >= extent;
            p_corr = 1 - (1 - p_uncorr) .^ sum(which);

            for i = 1:length(p_uncorr)
                if ~which(i)
                    continue;
                end
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
