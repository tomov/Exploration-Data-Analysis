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
            figure('pos', [10 10 660 700]);

            subplot(3,2,1);
            imshow('images/figure1Aleft.png'); %, 'InitialMagnification', 'fit');  

            fontsize = 16;
            linewidth = 3;
            markersize = 6;

            %% Sam_Figure1
            subplot(3,2,2);

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

            %% Sam_Figure2

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
            
            T = {'Relative uncertainty: intercept shift' 'Total uncertainty: slope shift'};
            for i = 1:2
                subplot(3,2,i+2);
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
                subplot(3,2,i+4);
                
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
            
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.06, 0.95, 'A', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.67, 'B', 'FontSize', 25, 'FontWeight', 'bold');
            text(0.06, 0.35, 'C', 'FontSize', 25, 'FontWeight', 'bold');
            

        case 'Figure2'
            % RU - trial contrast 
            %
            figure('pos', [100 100 0.75*653 0.75*352]);

            fontsize = 12;
           
            %PICpng = imread('images/badre_ROI.png');   %  <-- to cross-check 
            PICpng = imread('images/RU-trial.png');

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
            plot(pline_x, pline_y, '-', 'LineWidth', 2, 'Color', [1 1 1]);
            hold off;

            title('RU - trial', 'FontSize', fontsize);


        case 'Figure3'
            % TU - trial contrast
            %

            load cross_subject_glm21_TU_TU_-_trial_sphere.mat;

            figure('pos', [100 100 1.5*653 1.5*552]);


            fontsize = 12;
           
            h = subplot(1,2,1);
            imshow('images/TU-trial.png'); %, 'InitialMagnification', 'fit');  
            title('TU - trial', 'FontSize', fontsize);

            h = subplot(3,2,4);
            pos = get(h, 'position');
            pos(3) = pos(3) * 0.6;
            pos(1) = pos(1) * 0.9;
            subplot(3,2, 4, 'position', pos);
            scatter(all_b{6}', w(:,3));
            lsline;
            xlabel('\beta_{TU}');
            ylabel('w_3');
            title('Insula (L)', 'interpreter', 'none', 'FontSize', fontsize);

            str = sprintf('r = %.2f, p = %.3f', r(6), p_uncorr(6));
            text(-0.15,0.004, str, 'FontSize', 9);

            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            axes(ax1);
            text(0.13, 0.65, 'A', 'FontSize', 20, 'FontWeight', 'bold');
            text(0.48, 0.65, 'B', 'FontSize', 20, 'FontWeight', 'bold');


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
