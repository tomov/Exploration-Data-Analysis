function data = load_data(filename)
    
    % Load data.
    %
    % USAGE: data = load_data
    %
    % cond: 1 = RS, 2 = SR, 3 = RR, 4 = SS

    if ~exist('filename', 'var')
        filename = 'data.csv';
    end
    
    X = csvread(filename,1);
    %F = {'subject' 'run', 'block', 'trial', 'mu1', 'mu2', 'r1', 'r2', 'choice', 'reward', 'RT', 'cond', 'trial_onset', 'choice_onset', 'feedback_onset'};
    F = {'subject' 'block' 'trial' 'mu1' 'mu2' 'choice' 'reward' 'RT' 'cond'};

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
