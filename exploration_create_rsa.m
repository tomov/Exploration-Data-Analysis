function rsa = exploration_create_rsa(rsa_idx, subj, run)

    % Create rsa structure, helper function for creating EXPT in
    % exploration_expt.m
    %
    % USAGE: rsa = exploration_create_rsa(model,subj,run)
    %
    % INPUTS:
    %   rsa_idx - positive integer indicating which RSA we're doing
    %   subj - integer specifying which subject is being analyzed
    %   run - integer specifying the run
    %
    % OUTPUTS:
    %   rsa - a structure with the following fields:
    %     .event 
    %     .betas_glmodel
    %     .radius
    %     .model
    %
    % Momchil Tomov, Sep 2018


    fprintf('rsa %d, subj %d, run %d\n', rsa_idx, subj, run);

    data = load_data;
    conds = {'RS', 'SR', 'RR', 'SS'};
    

    [allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();
    
  
    % skip bad runs
    runs = find(goodRuns{subj});
    run = runs(run);
    fprintf('run %d \n', run);
    
    which_trials = data(subj).run == run;  
   
    fprintf('which_trials = %s\n', sprintf('%d', which_trials));

    % RSAs
    %
    switch rsa_idx

        % basic 
        %
        case 1
            rsa.event = 'trial_onset';
            rsa.betas_glmodel = 144;
            rsa.radius = 2.6666;

            rsa.model(1).name = 'condition';
            rsa.model(1).features = data(subj).cond(which_trials);
            rsa.model(1).distance_measure = @(c1, c2) c1 ~= c2;
            rsa.model(1).is_control = false;

            rsa.model(2).name = 'RS_SR_vs_RR_SS';
            rsa.model(2).features = data(subj).cond(which_trials) == 1 | data(subj).cond(which_trials) == 2;
            rsa.model(2).distance_measure = @(c1, c2) c1 ~= c2;
            rsa.model(2).is_control = false;

        % Sam's model
        %
        case 2
            [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'left');

            rsa.event = 'trial_onset';
            rsa.betas_glmodel = 144;
            rsa.radius = 2.6666;

            rsa.model(1).name = 'Qs_and_sigmas';
            rsa.model(1).features = [Q1, Q2, std1, std2];
            rsa.model(1).distance_measure = 'cosine';
            rsa.model(1).is_control = false;

            % controls

            rsa.model(2).name = 'time';
            rsa.model(2).features = data(subj).trial_onset(which_trials);
            rsa.model(2).distance_measure = 'euclidean';
            rsa.model(2).is_control = true;

            rsa.model(3).name = 'run';
            rsa.model(3).features = data(subj).run(which_trials);
            rsa.model(3).distance_measure = @(c1, c2) c1 ~= c2;
            rsa.model(3).is_control = true;

        otherwise
            assert(false, 'invalid rsa_idx -- should be one of the above');

    end % end of switch statement

end
