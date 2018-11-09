function show_figure(fig)

    switch fig
            
        case 'badre_ROI'
            % show the Badre ROI 
            %
            EXPT = exploration_expt();
            struc = fullfile(EXPT.modeldir,'mean.nii');
            masks = badre_2012_create_masks(false);
            bspmview(masks{1}, struc);


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
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference (V)','FontSize',fontsize);
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
                ylabel('Choice probability','FontSize',fontsize);
                xlabel('Expected value difference','FontSize',fontsize);
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
            
            % elot results
            [beta,~,stats] = fixedEffects(results);
            subplot(1,3,1);
            %errorbar(beta(1:4),stats.SE(1:4),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            errorbar(beta(1:4),(stats.Upper(1:4) - stats.Lower(1:4))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5],'YLim', [-1 1]);
            ylabel('Intercept','FontSize',fontsize);
            ylim([-0.5 0.5]);

            subplot(1,3,2);
            %errorbar(beta(5:8),stats.SE(5:8),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            errorbar(beta(5:8),(stats.Upper(5:8) - stats.Lower(5:8))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            set(gca,'FontSize',fontsize,'XTickLabel',{'RS' 'SR' 'RR' 'SS'},'XLim',[0.5 4.5], 'YLim', [0 3]);
            ylabel('Slope','FontSize',fontsize);
            ylim([1 2.5]);


            subplot(1,3,3);

            % Probit analysis of computational variables
            load results_glme_fig3
            results = results_VTURU;
                        
            % plot results
            [beta,~,stats] = fixedEffects(results);
            %errorbar(beta([3 1 2]),stats.SE([3 1 2]),'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            errorbar(beta([3 1 2]),(stats.Upper([3 1 2]) - stats.Lower([3 1 2]))/2,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
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

            load('main_effect_glm21_TU_RU_-_trial.mat');
            beta(2) = m(2);
            ci(2) = (cis{2}(2) - cis{2}(1)) / 2;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

          
            subplot(1,4,4);

            %{
            load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_RU_RU_-_trial_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{2});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            load('univariate_decoder_glm21_TU_RU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_TU_RU_-_trial_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{2});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Figure3C.mat', 'beta', 'err', 'ci');
            %}

            load Figure3C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\hat{RU}$', '$V/\hat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);
            


            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Figure3', '-dpdf');








        case 'Figure4'
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
            err(1) = (cis{end}(2) - cis{end}(1)) / 2;

            load('main_effect_glm21_TU_TU_-_trial.mat');
            beta(2) = m(end);
            err(2) = (cis{end}(2) - cis{end}(1)) / 2;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('Insula (L)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

         
            
            subplot(1,4,4);

            %{
            load('univariate_decoder_glm21_RU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_RU_TU_-_trial_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{end});
            beta(1) = b(4);
            err(1) = s.SE(4);
            ci(1) = (s.Upper(4) - s.Lower(4)) / 2;

            load('univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000_standardize=2_mixed=0.mat');
            %load('univariate_decoder_glm21_TU_TU_-_trial_norm=4_orth=1_lambda=1.000000.mat');
            [b,~,s] = fixedEffects(results_both{end});
            beta(2) = b(4);
            err(2) = s.SE(4);
            ci(2) = (s.Upper(4) - s.Lower(4)) / 2;

            save('Figure4C.mat', 'beta', 'err', 'ci');
            %}

            load Figure4C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('Insula (L)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\hat{RU}$', '$V/\hat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-11 4]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Figure4', '-dpdf');




        case 'Supp_Figure1'
            % Badre RLPFC 
            % parallels Figure 3
            %
            
            % EXPT = exploration_expt()
            % struc = fullfile(EXPT.modeldir, 'mean.nii')
            % masks = badre_2012_create_masks(false);
            % bspmview(masks{1}, struc);
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


            subplot(1,4,3);

            load('main_effect_glm21_RU_badre.mat');
            beta(1) = m(1);
            ci(1) = (cis{1}(2) - cis{1}(1)) / 2;

            load('main_effect_glm21_TU_badre.mat');
            beta(2) = m(1);
            ci(2) = (cis{1}(2) - cis{1}(1)) / 2;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

          
            subplot(1,4,4);

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

            load Supp_Figure1C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('RLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\hat{RU}$', '$V/\hat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 16]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);
            

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');

            print('images/Supp_Figure1', '-dpdf');


            


        case 'Supp_Figure2'
            % Badre DLPFC
            % Parallels Figure 4
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


            subplot(1,4,3);

            load('main_effect_glm21_RU_dlpfc.mat');
            beta(1) = m(1);
            err(1) = (cis{1}(2) - cis{1}(1)) / 2;

            load('main_effect_glm21_TU_dlpfc.mat');
            beta(2) = m(1);
            err(2) = (cis{1}(2) - cis{1}(1)) / 2;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,err,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('DLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2], 'XTickLabel',{'$|RU|$', '$TU$'},'XLim',[0.5 2.5], 'Ylim', [-0.1 0.2]);
            ylabel('Neural coefficient (\beta)','FontSize',fontsize);

         
            
            subplot(1,4,4);

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

            load Supp_Figure2C;

            plot([0 3],[0 0],'--','LineWidth',linewidth,'Color',[0.6 0.6 0.6]);
            hold on;
            errorbar(beta,ci,'ok','MarkerSize',markersize,'MarkerFaceColor','k');
            hold off;
            title('DLPFC (R)');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fontsize,'XTick', [1 2],'XTickLabel',{'$\hat{RU}$', '$V/\hat{TU}$'},'XLim',[0.5 2.5], 'Ylim', [-3 6]);
            ylabel('Regression coefficient (w)','FontSize',fontsize);

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.10, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.70, 0.95, 'C', 'FontSize', 20, 'FontWeight', 'bold');


            print('images/Supp_Figure2', '-dpdf');





        case 'Supp_Figure3'
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


            print('images/Supp_Figure5', '-dpdf');



        case 'Supp_Table1'
            % RU - trial

            ccnl_view(exploration_expt(), 21, 'RU - trial');
            ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'RU - trial', 0.001, '+/-', 0.05, 20, 1);

        case 'Table2'
            % TU - trial with p's

            ccnl_view(exploration_expt(), 21, 'TU - trial');
            ccnl_results_table('AAL2', 'peak', exploration_expt(), 21, 'TU - trial', 0.001, '+/-', 0.05, 20, 1);

           % load ...












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
