function data = load_data
    
    % Load data.
    %
    % USAGE: data = load_data
    %
    % cond: 1 = RS, 2 = SR, 3 = RR, 4 = SS
    
    X = csvread('data.csv',1);
    F = {'subject' 'run', 'block', 'trial', 'mu1', 'mu2', 'r1', 'r2', 'choice', 'reward', 'RT', 'cond', 'trial_onset', 'choice_onset', 'feedback_onset'};
    S = unique(X(:,1));
    
    for s = 1:length(S)
        ix = X(:,1)==S(s);
        for f = 1:length(F)
            data(s).(F{f}) = X(ix,f);
        end
        
        [~,k] = max([data(s).mu1 data(s).mu2],[],2);
        data(s).acc = mean(data(s).choice==k);
        data(s).timeout = isnan(data(s).choice);
    end
    
    acc = [data.acc];
    data(acc<0.55) = [];

    % take care of timeouts
    % TODO do it properly
    %
    for s = 1:length(S)
        data(s).choice(data(s).timeout) = 1 + (rand(size(data(s).choice(data(s).timeout))) > 0.5); % random choices
        data(s).RT(data(s).timeout) = 2; % = choiceDuration = timeout
        data(s).choice_onset(data(s).timeout) = data(s).trial_onset(data(s).timeout) + data(s).RT(data(s).timeout);
    end
