function tbl = data2table(data,standardize,no_timeouts)
    
    if nargin < 2; standardize = 0; end
    if nargin < 3; no_timeouts = 1; end

    % take care of timeouts
    %
    if no_timeouts
        % get rid of them
        F = fieldnames(data);
        for s = 1:length(data)
            which_trials = ~data(s).timeout;
            for i = 1:length(F)
                if size(data(s).(F{i}), 1) == length(which_trials)
                    data(s).(F{i}) = data(s).(F{i})(which_trials,:);
                end
            end
        end
    else
        % assign random choices
        for s = 1:length(data)
            data(s).choice(data(s).timeout) = 1 + (rand(size(data(s).choice(data(s).timeout))) > 0.5); % random choices;
            data(s).RT(data(s).timeout) = 2; % = choiceDuration = timeout
            data(s).choice_onset(data(s).timeout) = data(s).trial_onset(data(s).timeout) + data(s).RT(data(s).timeout);
        end
    end
    
    RS = []; SS = []; C = []; V = []; S = []; RU = []; TU = []; rt = []; risky = []; cond = [];
    for s = 1:length(data)
        latents = kalman_filter(data(s));
        V = [V; latents.m(:,1) - latents.m(:,2)];
        risky = [risky; double((data(s).cond==1&data(s).choice==1)|(data(s).cond==2&data(s).choice==2)) - double((data(s).cond==1&data(s).choice==2)|(data(s).cond==2&data(s).choice==1))];
        RS = [RS; double(data(s).cond==1) - double(data(s).cond==2)];
        SS = [SS; double(data(s).cond==4) - double(data(s).cond==3)];
        cond = [cond; data(s).cond];
        C = [C; double(data(s).choice==1)];
        N = length(data(s).choice);
        S = [S; zeros(N,1)+s];
        TU = [TU; sqrt(latents.s(:,1) + latents.s(:,2))];
        RU = [RU; sqrt(latents.s(:,1)) - sqrt(latents.s(:,2))];
        rt = [rt; data(s).RT];
    end
    
    SSV = SS.*V;
    VTU = V./TU;
    rt = log(rt);
    cond = categorical(cond);
    
    if standardize==1
        VTU = zscore(VTU);
        V = zscore(V);
        RU = zscore(RU);
    end
    
    tbl = table(RS,SS,SSV,C,S,RU,VTU,V,rt,risky,cond);
