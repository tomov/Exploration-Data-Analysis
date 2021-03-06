% WRONG -- this is super circular;
% the RU and TU ROIs correlate with RU and TU
% the DV ROI correlates with DV which correlates with RU and TU
% so ofc their beta series will be correlated doh

% correlate beta series between DV and RU and TU ROIs
%
function beta_series_corr()

    printcode;

    DV_glm = 47;
    VTURU_glm = 36;
    beta_series_glm = 23;


    data = load_data;

    EXPT = exploration_expt();
    clusterFWEcorrect = false;
    extent = 100;
    Num = 1;

    load results_glme_fig3_nozscore.mat;
    w = getEffects(results_VTURU, false);

    % get DV ROI from GLM 47
    [DV_masks, ~] = get_masks(DV_glm, 'DV', clusterFWEcorrect, extent, Num);

    % get RU ROIs from GLM 36
    [RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);

    % get TU ROIs from GLM 36
    [TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);


    seed_masks = [RU_masks, TU_masks];


    for c = 1:length(seed_masks)
        seed_mask = seed_masks{c};

        fprintf('\n\n seed mask = %s \n\n', seed_mask);

        % RU beta series
        for s = 1:length(data)
            betas{s} = get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', seed_mask);
        end

        clear R;

        mean_rs = [];
        ttest_ts = [];
        ttest_ps = [];
        for dv = 1:length(DV_masks)

            DV_mask = DV_masks{dv};

            % DV beta series
            for s = 1:length(data)
                DV_betas{s} = get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', DV_mask);
            end

            rs = [];
            for s = 1:length(data)
                [r, p] = corr(betas{s}, DV_betas{s});
                rs(s) = r;
            end

            rs
            mean_rs(dv,:) = mean(rs);

            % is coupling > 0 across subjects?
            %
            rs = atanh(rs);
            [h, p, ci, stats] = ttest(rs)

            ttest_ts(dv,:) = stats.tstat;
            ttest_ps(dv,:) = p;


            % correlate f'n coupling with w's
            %
            if c == 1
                % TODO hardcoded
                [r, p] = corr(rs', w(:,2));
            else
                [r, p] = corr(rs', w(:,3));
            end

            cc_rs(dv,:) = r;
            cc_ps(dv,:) = p;
        end

        T = table(DV_masks', mean_rs, ttest_ts, ttest_ps, cc_rs, cc_ps);
        disp(T);
        all_T{c} = T;
    end

    save('beta_series_corr.mat', '-v7.3');
end
