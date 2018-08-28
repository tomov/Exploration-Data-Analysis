% helper f'n wrapper around kalman_filter
%
% how = how to pick option 1 so as to disambiguate option 1 vs. option 2 and determine the sign of RU, V, etc.
%     how = left -> 1 = left, 2 = right
%     how = max -> 1 = max, 2 = min
%     how = chosen -> 1 = chosen, 2 = unchosen
%     how = abs -> return abs value; option varies depending on quantity
%
function [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, how)
    latents = kalman_filter(data(subj));
   
    %load results_glme_fig3_nozscore_withtimeouts.mat; <--- the old one :( wrong
    load results_glme_fig3_nozscore.mat;
 
    [w_f, names_f] = fixedEffects(results_VTURU);
    [w_r, names_r] = randomEffects(results_VTURU);

    w = w_f + w_r((subj - 1) * 3 + 1 : subj * 3);

    w1 = w(3);
    w2 = w(1);
    w3 = w(2);
    w = [w1 w2 w3]; % reorder them

    QL = latents.m(which_trials,1);
    QR = latents.m(which_trials,2);
    stdL = sqrt(latents.s(which_trials,1));
    stdR = sqrt(latents.s(which_trials,2));

    % sanity check
    V = QL - QR;
    RU = stdL - stdR;
    TU = sqrt(stdL.^2 + stdR.^2);
    DV = w1 * V + w2 * RU + w3 * V./TU;
    pred = normcdf(DV); % manual prediction
    tbl = data2table(data,0,0); % don't exclude timeouts here b/c we're predicting
    tbl = tbl(table2array(tbl(:,'S')) == subj, :);
    y = predict(results_VTURU, tbl);
    y = y(which_trials); % glm prediction
    assert(immse(y, pred) < 1e-10, 'GLM prediction different from computed prediction');

    
    TU = sqrt(stdL.^2 + stdR.^2);
    RU = nan(size(TU));
    V = nan(size(TU));
    Q1 = nan(size(TU));
    std1 = nan(size(TU));
    Q2 = nan(size(TU));
    std2 = nan(size(TU));
    DV = nan(size(TU)); % decision value
    DQ1 = nan(size(TU)); % decision value for option 1
    DQ2 = nan(size(TU)); % decision value for option 2
    
    DQL = (w1*QL) + (w2*stdL) + ((QL./TU)*w3);  
    DQR = (w1*QR) + (w2*stdR) + ((QR./TU)*w3);
    assert(immse(y, normcdf(DQL - DQR)) < 1e-10, 'GLM prediction different from computed prediction');

    switch how

        case 'left' % 1 = left, 2 = right

           DQ1 = DQL;
           DQ2 = DQR;
 
           Q1 = QL;
           Q2 = QR;
 
           std1 = stdL;
           std2 = stdR;

        case 'max' % 1 = max, 2 = min

            for i=1:length(DQL)
               if DQL(i) >= DQR(i)
                   DQ1(i) = DQL(i);
                   DQ2(i) = DQR(i);
         
                   Q1(i) = QL(i);
                   Q2(i) = QR(i);
         
                   std1(i) = stdL(i);
                   std2(i) = stdR(i);
               else 
                   DQ1(i) = DQR(i);
                   DQ2(i) = DQL(i);
         
                   Q1(i) = QR(i);
                   Q2(i) = QL(i);
         
                   std1(i) = stdR(i);
                   std2(i) = stdL(i);
               end 
            end

        case 'chosen' % 1 = chosen, 2 = unchosen
            assert(~any(which_trials & data(subj).timeout), 'Do not include timeouts when how = chosen');

            for i=1:length(DQL)
               if data(subj).choice(i) == 1 % chose left
                   DQ1(i) = DQL(i);
                   DQ2(i) = DQR(i);
         
                   Q1(i) = QL(i);
                   Q2(i) = QR(i);
         
                   std1(i) = stdL(i);
                   std2(i) = stdR(i);
               else 
                   DQ1(i) = DQR(i);
                   DQ2(i) = DQL(i);
         
                   Q1(i) = QR(i);
                   Q2(i) = QL(i);
         
                   std1(i) = stdR(i);
                   std2(i) = stdL(i);
               end 
           end

        case 'abs' % 1 = bigger, 2 = smaller

            for i=1:length(DQL)
                DQ1(i) = max(DQL(i), DQR(i));
                DQ2(i) = min(DQL(i), DQR(i));
     
                Q1(i) = max(QL(i), QR(i));
                Q2(i) = min(QL(i), QR(i));
     
                std1(i) = max(stdL(i), stdR(i));
                std2(i) = min(stdL(i), stdR(i));
            end

        otherwise
            assert(false, 'Invalid parameter for how');
    end
   
    RU = std1 - std2; 
    V = Q1 - Q2;
    DV = DQ1 - DQ2;
         
   
    %getMax = max(latents.s(:,1),latents.s(:,2));
    %getMin = max(latents.s(:,1),latents.s(:,2));
    
    
    %RU = sqrt(latents.s(:,1)) - sqrt(latents.s(:,2));
    %V = latents.m(:,1) - latents.m(:,2);
    VTU = V./TU;
    assert(all(size(VTU) == size(QL)));
end
