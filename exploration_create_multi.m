function multi = exploration_create_multi(glmodel, subj, run, save_output)

    % Create multi structure, helper function for creating EXPT in
    % exploration_expt.m
    %
    % USAGE: multi = exploration_create_multi(model,subj,run)
    %
    % INPUTS:
    %   glmodel - positive integer indicating general linear model
    %   subj - integer specifying which subject is being analyzed
    %   run - integer specifying the run
    %

    % OUTPUTS:
    %   multi - a sctructure with the folloowing fields
    %        .names{i}
    %        .onsets{i}
    %        .duratlsions{i}
    %        optional:
    %        .pmod(i).name
    %        .pmod(i).param
    %        .pmod(i).poly
    %
    % Momchil Tomov, July 2018

    if nargin < 4 || isempty(save_output)
        save_output = false;
    end

    fprintf('glm %d, subj %d, run %d\n', glmodel, subj, run);

    data = load_data;

    % handle timeouts -- anecdotally, people pressed even if there was a timeout (just too late) -> count press at timeout
    for s = 1:length(data)
        data(s).RT(data(s).timeout) = 2; % = choiceDuration = timeout
        data(s).choice_onset(data(s).timeout) = data(s).trial_onset(data(s).timeout) + data(s).RT(data(s).timeout);
    end
    
    conds = {'RS', 'SR', 'RR', 'SS'};
    
    

    [allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();
   
    SPM_run = run; % save the SPM run (SPM doesn't see bad runs)
  
    % skip bad runs
    runs = find(goodRuns{subj});
    run = runs(run);
    fprintf('run %d \n', run);
    

    % distribution = results_V.Distribution;
    % dispersion = results_V.Dispersion;
    % dispersion_estimate = results_V.DispersionEstimate;
    % num_variables = results_V.NumVariables;
    % num_predictors = results_V.NumPredictors;
    % num_observation = results_V.NumObservation;
    % variables = results_V.Variables;
    % loglik    = results_V.LogLikelihood;
    % dfe = results_V.DFE;
    % sse = results_V.SSE;
    % sst = results_V.SST;
    % ssr = results_V.SSR;
    % coefcov = results_V.CoefficientCovariance;
    
    which_trials = data(subj).run == run;  
   
    fprintf('which_trials = %s\n', sprintf('%d', which_trials));


    %[results_V, results_VTU, results_VRU, results_VTURU ] = model_comparison(data(subj));


    % GLMs
    %
    switch glmodel

        % #Sam
        % Condition @ trial_onset 
        % nuisance @ choice and feedback onset 
        % 
        case 1 
            % condition @ trial onset
            %
            for cond = 1:4
                multi.names{cond} = conds{cond};
                multi.onsets{cond} = data(subj).trial_onset(which_trials & data(subj).cond == cond)';
                multi.durations{cond} = zeros(size(multi.onsets{cond}));
            end

            % nuisance @ choice onset
            %
            multi.names{5} = 'choice_onset';
            multi.onsets{5} = data(subj).choice_onset(which_trials)';
            multi.durations{5} = zeros(size(multi.onsets{5}));

            % nuisance @ feedback onset
            %
            multi.names{6} = 'feedback_onset';
            multi.onsets{6} = data(subj).feedback_onset(which_trials)';
            multi.durations{6} = zeros(size(multi.onsets{6}));
         
         
        % #Sam
        % glm 1 with longer durations
        % Condition @ trial_onset, duration = RT
        % nuisance @ choice and feedback onset 
        % 
        case 2
            % condition @ trial onset
            %
            for cond = 1:4
                which = which_trials & data(subj).cond == cond;

                multi.names{cond} = conds{cond};
                multi.onsets{cond} = data(subj).trial_onset(which)';
                multi.durations{cond} = data(subj).RT(which)' / 1000;
            end

            % nuisance @ choice onset
            %
            multi.names{5} = 'choice_onset';
            multi.onsets{5} = data(subj).choice_onset(which_trials)';
            multi.durations{5} = zeros(size(multi.onsets{5}));

            % nuisance @ feedback onset
            %
            multi.names{6} = 'feedback_onset';
            multi.onsets{6} = data(subj).feedback_onset(which_trials)';
            multi.durations{6} = zeros(size(multi.onsets{6}));
       

        % #Sam
        % glm 1 with longer durations
        % Condition @ trial_onset, duration = trial duration
        % nuisance @ choice and feedback onset 
        % 
        case 3
            % condition @ trial onset
            %
            for cond = 1:4
                which = which_trials & data(subj).cond == cond;

                multi.names{cond} = conds{cond};
                multi.onsets{cond} = data(subj).trial_onset(which)';
                multi.durations{cond} = (1 + data(subj).feedback_onset(which) - data(subj).trial_onset(which))';
            end

            % nuisance @ choice onset
            %
            multi.names{5} = 'choice_onset';
            multi.onsets{5} = data(subj).choice_onset(which_trials)';
            multi.durations{5} = zeros(size(multi.onsets{5}));

            % nuisance @ feedback onset
            %
            multi.names{6} = 'feedback_onset';
            multi.onsets{6} = data(subj).feedback_onset(which_trials)';
            multi.durations{6} = zeros(size(multi.onsets{6}));
        


        % ---------------------- max - min ---------------------------------------x


        % #Sam
        % RU, TU @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 4
           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
       

        % #Sam
        % RU, TU, V, V/TU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 5
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
        
               
        % #Sam
        % decision value @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 6
           [~, ~, ~, ~, DV] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % -----------------------------   chosen - unchosen ------------------------------------------------------------


        % #Sam
        % glm 4 but chosen - unchosen
        % RU, TU @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 7
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
       

        % #Sam
        % glm 5 but chosen - unchosen
        % RU, TU, V, V/TU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 8
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
        
               
        % #Sam
        % glm 6 but chosen - unchosen
        % decision value @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 9
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [~, ~, ~, ~, DV] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % ------------------------------ |abs values| --------------------------------------


        % #Sam
        % glm 4 but with |abs|
        % RU, TU @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 10
           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
       

        % #Sam
        % glm 5 but with |abs|
        % RU, TU, V, V/TU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        % REVIEWER #1 -- odd runs only! see exploration_getSubjectsDirsAndRuns
        % for ROI selection
        %
        case 11
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           assert(mod(run, 2) == 1); % make sure we didn't fuck up in exploration_getSubjectsDirsAndRuns

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
        
               
        % #Sam ; IGNORE -- identical to glm 6
        % glm 6 but with |abs|
        % decision value @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 12
           [~, ~, ~, ~, DV] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % ---------------------------------------------------------------------------------


        % #Sam
        % same as 4 but with trial # pmod
        % RU, TU, trial # @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 13
           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % #Sam 
        % TU, split DV and RU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 14
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'std1';
           multi.pmod(1).param{1} = std1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'std2';
           multi.pmod(1).param{2} = std2';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'DQ1';
           multi.pmod(1).param{3} = DQ1';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'DQ2';
           multi.pmod(1).param{4} = DQ2';
           multi.pmod(1).poly{4} = 1; 

           multi.pmod(1).name{5} = 'TU';
           multi.pmod(1).param{5} = TU';
           multi.pmod(1).poly{5} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % Q1+std1, Q2+std2 @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 15
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'max');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'Q1std1';
           multi.pmod(1).param{1} = w(1) * Q1' + w(2) * std1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'Q2std2';
           multi.pmod(1).param{2} = w(1) * Q2' + w(2) * std2';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % glm 15 but chosen - unchosen
        % Q1+std1, Q2+std2 @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 16
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'Q1std1';
           multi.pmod(1).param{1} = w(1) * Q1' + w(2) * std1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'Q2std2';
           multi.pmod(1).param{2} = w(1) * Q2' + w(2) * std2';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % #Sam 
        % TU, split DV and RU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 17
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'std1';
           multi.pmod(1).param{1} = std1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'std2';
           multi.pmod(1).param{2} = std2';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'DQ1';
           multi.pmod(1).param{3} = DQ1';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'DQ2';
           multi.pmod(1).param{4} = DQ2';
           multi.pmod(1).poly{4} = 1; 

           multi.pmod(1).name{5} = 'TU';
           multi.pmod(1).param{5} = TU';
           multi.pmod(1).poly{5} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

 
        % glm 15 but left - right
        % Q1+std1, Q2+std2 @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 18
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'left');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'Q1std1';
           multi.pmod(1).param{1} = w(1) * Q1' + w(2) * std1';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'Q2std2';
           multi.pmod(1).param{2} = w(1) * Q2' + w(2) * std2';
           multi.pmod(1).poly{2} = 1;

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % #Sam
        % glm 10 with trial # pmod
        % |RU|, TU, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 19
           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));
      

        % glm 19 w/o TU
        % |RU|, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 20
           [~, RU, TU, ~] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'trial';
           multi.pmod(1).param{2} = data(subj).trial(which_trials)';
           multi.pmod(1).poly{2} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % #Sam
        % glm 19 without timeouts
        % |RU|, TU, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 21
           [~, RU, TU, ~] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials & ~data(subj).timeout)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % glm 18 but without trial_onset_L
        % Q1+std1, Q2+std2 @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 22
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'left');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'Q1std1';
           multi.pmod(1).param{1} = w(1) * Q1' + w(2) * std1';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'Q2std2';
           multi.pmod(1).param{2} = w(1) * Q2' + w(2) * std2';
           multi.pmod(1).poly{2} = 1;

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

        % for RSA
        % similar to GLM 143 in ConLearn
        %
        case 23
           idx = 0;

           block = data(subj).block(which_trials);
           trial = data(subj).trial(which_trials);

           onsets = data(subj).trial_onset(which_trials);
           for t = 1:numel(onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_block_', num2str(block(t)), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['trial_onset_', suffix];
               multi.onsets{idx} = [onsets(t)];
               multi.durations{idx} = [0];
           end

           onsets = data(subj).choice_onset(which_trials);
           for t = 1:numel(onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_block_', num2str(block(t)), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['choice_onset_', suffix];
               multi.onsets{idx} = [onsets(t)];
               multi.durations{idx} = [0];
           end

           onsets = data(subj).feedback_onset(which_trials);
           for t = 1:numel(onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_block_', num2str(block(t)), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['feedback_onset_', suffix];
               multi.onsets{idx} = [onsets(t)];
               multi.durations{idx} = [0];
           end

           
        % #Sam
        % glm 21 with RU = chosen - unchosen
        % RU, TU, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 24
           [~, RU, TU, ~] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'chosen'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials & ~data(subj).timeout)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % #Sam
        % glm 21 with RU = left - right 
        % RU, TU, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 25
           [~, RU, TU, ~] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'left'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials & ~data(subj).timeout)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % #Sam
        % glm 21 with RU = max - min
        % RU, TU, trial pmod @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 26
           [~, RU, TU, ~] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'max'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'trial';
           multi.pmod(1).param{3} = data(subj).trial(which_trials & ~data(subj).timeout)';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % #Sam
        % D1 + D2 = w1*Q1 + w2*std1 + w1*Q2 + w2*std2 @ trial onset
        % no timeouts
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        % nuisance @ trial onset for timeouts
        %
        case 27
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'left');
           D1 = w(1) * Q1 + w(2) * std1;
           D2 = w(1) * Q2 + w(2) * std2;
           both = D1 + D2;

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = both';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end

        % #Sam similar to 27
        % |D1 - D2| = |w1*Q1 + w2*std1 - w1*Q2 - w2*std2| @ trial onset
        % no timeouts
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        % nuisance @ trial onset for timeouts
        %
        case 28
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'left');
           D1 = w(1) * Q1 + w(2) * std1;
           D2 = w(1) * Q2 + w(2) * std2;
           diff = abs(D1 - D2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'diff';
           multi.pmod(1).param{1} = diff';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % #Sam
        % glm 12 without timeouts
        % |DV| @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        % goes w/ 45
        %
        case 29
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end

        % #Sam
        % V without timeouts
        % |V| @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 30
           [V, RU, TU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'V';
           multi.pmod(1).param{1} = V';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % Q's, std's, Q/TU's @ trial_onset (1 = left, 2 = right)
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 31
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'left');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'Q1';
           multi.pmod(1).param{1} = Q1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'Q2';
           multi.pmod(1).param{2} = Q2';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'std1';
           multi.pmod(1).param{3} = std1';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'std2';
           multi.pmod(1).param{4} = std2';
           multi.pmod(1).poly{4} = 1; 

           multi.pmod(1).name{5} = 'Q1TU';
           multi.pmod(1).param{5} = Q1' ./ TU';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'Q2TU';
           multi.pmod(1).param{6} = Q2' ./ TU';
           multi.pmod(1).poly{6} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % Q's, std's, Q/TU's @ trial_onset (1 = bigger, 2 = smaller)
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 32
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'Q1';
           multi.pmod(1).param{1} = Q1';
           multi.pmod(1).poly{1} = 1; 

           if subj == 17 && run == 3
               % colinear regressors...
               multi.pmod(1).name{2} = 'Q2';
               multi.pmod(1).param{2} = Q2' + rand(size(Q2')) * 0.01;
               multi.pmod(1).poly{2} = 1; 
           else 
               multi.pmod(1).name{2} = 'Q2';
               multi.pmod(1).param{2} = Q2';
               multi.pmod(1).poly{2} = 1; 
           end

           multi.pmod(1).name{3} = 'std1';
           multi.pmod(1).param{3} = std1';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'std2';
           multi.pmod(1).param{4} = std2';
           multi.pmod(1).poly{4} = 1; 

           multi.pmod(1).name{5} = 'Q1TU';
           multi.pmod(1).param{5} = Q1' ./ TU';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'Q2TU';
           multi.pmod(1).param{6} = Q2' ./ TU';
           multi.pmod(1).poly{6} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % Q's, std's, Q/TU's @ trial_onset (1 = chosen, 2 = unchosen)
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 33
           which_trials = which_trials & ~data(subj).timeout; % exclude timeouts
           fprintf('which_trials = %s\n', sprintf('%d', which_trials));

           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, 'chosen');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'Q1';
           multi.pmod(1).param{1} = Q1';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'Q2';
           multi.pmod(1).param{2} = Q2';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'std1';
           multi.pmod(1).param{3} = std1';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'std2';
           multi.pmod(1).param{4} = std2';
           multi.pmod(1).poly{4} = 1; 

           multi.pmod(1).name{5} = 'Q1TU';
           multi.pmod(1).param{5} = Q1' ./ TU';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'Q2TU';
           multi.pmod(1).param{6} = Q2' ./ TU';
           multi.pmod(1).poly{6} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % #Sam
        % same as 11 but by condition
        % see also 3
        %
        case 34
           % condition @ trial onset
           %
           for cond = 1:4
               which = which_trials & data(subj).cond == cond;
               [V, RU, TU, VTU] = get_latents(data, subj, which, 'abs');

               multi.names{cond} = conds{cond};
               multi.onsets{cond} = data(subj).trial_onset(which)';
               multi.durations{cond} = zeros(size(multi.onsets{cond}));

               multi.orth{cond} = 0; % do not orthogonalise them  

               multi.pmod(cond).name{1} = 'RU';
               multi.pmod(cond).param{1} = RU';
               multi.pmod(cond).poly{1} = 1;    

               multi.pmod(cond).name{2} = 'TU';
               multi.pmod(cond).param{2} = TU';
               multi.pmod(cond).poly{2} = 1; 

               % VTU and V introduce colinearity
           end

           multi.names{1 + cond} = 'choice_onset';
           multi.onsets{1 + cond} = data(subj).choice_onset(which_trials);
           multi.durations{1 + cond} = zeros(size(multi.onsets{1 + cond}));

           multi.names{2 + cond} = 'feedback_onset';
           multi.onsets{2 + cond} = data(subj).feedback_onset(which_trials);
           multi.durations{2 + cond} = zeros(size(multi.onsets{2 + cond}));

           multi.names{3 + cond} = 'trial_onset_L';
           multi.onsets{3 + cond} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{3 + cond} = zeros(size(multi.onsets{3 + cond}));


        % literally the same as GLM 11 but with even runs only (for further validating ROIs / testing)
        % REVIEWER #1 -- see exploration_getSubjectsDirsAndRuns
        %
        case 35
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           assert(mod(run, 2) == 0); % make sure we didn't fuck up in exploration_getSubjectsDirsAndRuns

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % #Sam
        % glm 5 but with |abs|
        % RU, TU, V, V/TU @ trial_onset 
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        % same as 11 and 35 but with all runs
        %
        case 36 % <----------------------------------------------------- THIS IS IT ------------------------------
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % exact copy of 36 but for ncf (Cent OS 6)
        %
        case 37
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % copy of 36 but across conditions
        %
        case 38

           for cond = 1:4
               which = which_trials & data(subj).cond == cond;

               [V, RU, TU, VTU] = get_latents(data, subj, which, 'abs');

               multi.names{cond} = conds{cond};
               multi.onsets{cond} = data(subj).trial_onset(which)';
               multi.durations{cond} = zeros(size(multi.onsets{cond}));

               multi.orth{cond} = 0; % do not orthogonalise them  

               multi.pmod(cond).name{1} = 'RU';
               multi.pmod(cond).param{1} = RU';
               multi.pmod(cond).poly{1} = 1;    

               multi.pmod(cond).name{2} = 'TU';
               multi.pmod(cond).param{2} = TU';
               multi.pmod(cond).poly{2} = 1; 

               multi.pmod(cond).name{3} = 'V';
               multi.pmod(cond).param{3} = V';
               multi.pmod(cond).poly{3} = 1; 

               multi.pmod(cond).name{4} = 'VTU';
               multi.pmod(cond).param{4} = VTU';
               multi.pmod(cond).poly{4} = 1; 

           end

           multi.names{5} = 'choice_onset';
           multi.onsets{5} = data(subj).choice_onset(which_trials);
           multi.durations{5} = zeros(size(multi.onsets{5}));

           multi.names{6} = 'feedback_onset';
           multi.onsets{6} = data(subj).feedback_onset(which_trials);
           multi.durations{6} = zeros(size(multi.onsets{6}));

           multi.names{7} = 'trial_onset_L';
           multi.onsets{7} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{7} = zeros(size(multi.onsets{7}));



        % copy of 36 but with RU only -- for VIFs followup 
        %
        case 39
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % copy of 36 but with TU only -- for VIFs followup 
        %
        case 40
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'TU';
           multi.pmod(1).param{1} = TU';
           multi.pmod(1).poly{1} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));




        % copy of 36 but with V only -- for VIFs followup 
        %
        case 41
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'V';
           multi.pmod(1).param{1} = V';
           multi.pmod(1).poly{1} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));




        % copy of 36 but with V/TU only -- for VIFs followup 
        %
        case 42
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'VTU';
           multi.pmod(1).param{1} = VTU';
           multi.pmod(1).poly{1} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % same as 29 but CentOS 7 (paired with 45)
        %
        % |DV| @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 43 % <-- nothing pos; L precentral negative 
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs'); % exclude timeouts

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % |V + RU| @ trial_onset ; proxy for DV
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 44 % <-- a bunch of negative blobs in visual and PPC
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'left'); % exclude timeouts

           DV = w(1) * V + w(2) * RU;
           DV = abs(DV);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); % exclude timeouts
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % same as 36 but without timeouts  UGH
        % goes w/ 29
        %
        case 45  %  <--- basically the same result; smaller ROIs for RU
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end

        % same as 36 but w/o choice_onset 
        %
        case 46 % <-- bigger R RLPFC (205), still n.s. after correction
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'feedback_onset';
           multi.onsets{2} = data(subj).feedback_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'trial_onset_L';
           multi.onsets{3} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{3} = zeros(size(multi.onsets{3}));



        % same as 43 but WITH timeouts, like 36
        %
        % |DV| @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 47  % <-- L precentral negative
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials); 
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % w1 * Q_chosen + w2 * std_chosen @ trial_onset; proxy for DV
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset 
        %
        case 48 % <-- parietal bilateral
           [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'chosen'); 

           DV = w(1) * Q1 + w(2) * std1;

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout); 
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % 36 but with residuals from DV Precentral (L) ROI (GLM 47)
        % for functional connectivity analysis (activation-induced correlations; see p. 133 from Poldrack book)
        %
        case 49 % WRONG -- see below; never ran it
           roi_glmodel = 47; % DV only
           roi_contrast = 'DV';
           clusterFWEcorrect = false;
           extent = 100;
           Num = 1;
           EXPT = exploration_expt();

           % get DV ROI from GLM 47
           [masks, region] = get_masks(roi_glmodel, roi_contrast, clusterFWEcorrect, extent, Num);
           DV_mask = masks{1};

           % extract residuals from DV ROI (from GLM 36!! we don't have DV there b/c of lineariy)
           res_glmodel = 36; % V, RU, TU, V/TU
           res = mean(ccnl_get_residuals(EXPT, res_glmodel, DV_mask, subj), 2);

           % subset residuals for this run only
           nTRs = 242;
           assert(mod(length(res), nTRs) == 0);
           res = res((SPM_run - 1)*nTRs + 1 : SPM_run*nTRs, :);
       
           % same as GLM 36
           %
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           % except add residuals from DV ROI !
           %
           assert(length(res) == nTRs);
           TR = EXPT.TR;
           % TODO must deconvolve with HRF... or add to raw design matrix, i.e. SPM.xX.X ... anyway, fuck this
           multi.names{2} = 'TR_onset';
           multi.onsets{2} = TR/2 : TR : nTRs*TR;
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.orth{2} = 0; % do not orthogonalise them  

           multi.pmod(2).name{1} = 'DVres';
           multi.pmod(2).param{1} = res';
           multi.pmod(2).poly{1} = 1;    

           % same nuisance regressors as 36
           %
           multi.names{3} = 'choice_onset';
           multi.onsets{3} = data(subj).choice_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'feedback_onset';
           multi.onsets{4} = data(subj).feedback_onset(which_trials);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           multi.names{5} = 'trial_onset_L';
           multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{5} = zeros(size(multi.onsets{5}));


        % same as 47 but at choice onset 
        % trying to make the activations positive
        %
        % |DV| @ choice_onset 
        % left choice @ choice_onset (!!!)
        % nuisance @ trial_onset and feedback_onset 
        %
        case 50  % <-- same; negatives
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'choice_onset';
           multi.onsets{1} = data(subj).choice_onset(which_trials); 
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'trial_onset';
           multi.onsets{2} = data(subj).trial_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'choice_onset_L';
           multi.onsets{4} = data(subj).choice_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % copy of 36, but left - right
        % seeking RU after FWE
        %
        case 51 % <-- RU - weird occipital, no RLPFC; TU same; nothing for V; absolutely nothing for VTU
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'left');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % 36 + beta series for RU and TU i.e. functional connectivity, WITH orth
        % => who are the RU and TU ROI's talking to?
        %
        case 52  % <-- RU_betas, TU_betas1 = lol whole brain ; TU_betas2 = frontoparietal 
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 1; % DO orthogonalise them -- need to be ultra-conservative here

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           %
           % brain activity pmods
           %

           VTURU_glm = 36;
           EXPT = exploration_expt();
           clusterFWEcorrect = false;
           extent = 100;
           Num = 1;

           beta_series_glm = 23;
           event = ['trial_onset_run_', num2str(run)];

           % get RU ROIs from GLM 36 and beta series from GLM 23
           RU_rois = [1];
           [RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);
           RU_mask = RU_masks{RU_rois(1)};

           RU_betas = get_beta_series(EXPT, beta_series_glm, subj, event, RU_mask);

           % get TU ROIs from GLM 36 and beta series from GLM 23
           TU_rois = [2, 8];
           [TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);
           TU_mask1 = TU_masks{TU_rois(1)};
           TU_mask2 = TU_masks{TU_rois(2)};

           TU_betas1 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask1);
           TU_betas2 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask2);

           multi.pmod(1).name{5} = 'RU_betas';
           multi.pmod(1).param{5} = RU_betas';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'TU_betas1';
           multi.pmod(1).param{6} = TU_betas1';
           multi.pmod(1).poly{6} = 1; 

           multi.pmod(1).name{7} = 'TU_betas2';
           multi.pmod(1).param{7} = TU_betas2';
           multi.pmod(1).poly{7} = 1; 


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));



        % same as 52 but do NOT orthogonalise
        %
        case 53  % <-- RU_betas = frontoparietal, TU_betas1 = whole brain ; TU_betas2 = frontoparietal ; RU = nothing positive => outcompeted by RU_betas => stick w/ 52
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do NOT orthogonalise them -- need to be ultra-conservative here

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           %
           % brain activity pmods
           %

           VTURU_glm = 36;
           EXPT = exploration_expt();
           clusterFWEcorrect = false;
           extent = 100;
           Num = 1;

           beta_series_glm = 23;
           event = ['trial_onset_run_', num2str(run)];

           % get RU ROIs from GLM 36 and beta series from GLM 23
           RU_rois = [1];
           [RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);
           RU_mask = RU_masks{RU_rois(1)};

           RU_betas = get_beta_series(EXPT, beta_series_glm, subj, event, RU_mask);

           % get TU ROIs from GLM 36 and beta series from GLM 23
           TU_rois = [2, 8];
           [TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);
           TU_mask1 = TU_masks{TU_rois(1)};
           TU_mask2 = TU_masks{TU_rois(2)};

           TU_betas1 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask1);
           TU_betas2 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask2);

           multi.pmod(1).name{5} = 'RU_betas';
           multi.pmod(1).param{5} = RU_betas';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'TU_betas1';
           multi.pmod(1).param{6} = TU_betas1';
           multi.pmod(1).poly{6} = 1; 

           multi.pmod(1).name{7} = 'TU_betas2';
           multi.pmod(1).param{7} = TU_betas2';
           multi.pmod(1).poly{7} = 1; 


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % 36 + condition regressor at trial onset
        % trying to get rid of negative activations for RU, also make it significant => could be that SS condition is screwing us (ppl stop paying attn)
        %
        case 54 % <-- even less stuff than 36 => less (regressors) is more (blobs)
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           for cond = 1:3 % skip last one => linearly dependent w/ trial_onset
                multi.names{cond + 4} = conds{cond};
                multi.onsets{cond + 4} = data(subj).trial_onset(which_trials & data(subj).cond == cond)';
                multi.durations{cond + 4} = zeros(size(multi.onsets{cond + 4}));
           end


        % 53 but with RU_betas only -- sanity check
        %
        case 55  % <--- whole brain; phew
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do NOT orthogonalise them -- need to be ultra-conservative here

           %
           % brain activity pmods
           %

           VTURU_glm = 36;
           EXPT = exploration_expt();
           clusterFWEcorrect = false;
           extent = 100;
           Num = 1;

           beta_series_glm = 23;
           event = ['trial_onset_run_', num2str(run)];

           % get RU ROIs from GLM 36 and beta series from GLM 23
           RU_rois = [1];
           [RU_masks, ~] = get_masks(VTURU_glm, 'RU', clusterFWEcorrect, extent, Num);
           RU_mask = RU_masks{RU_rois(1)};

           RU_betas = get_beta_series(EXPT, beta_series_glm, subj, event, RU_mask);

           % get TU ROIs from GLM 36 and beta series from GLM 23
           TU_rois = [2, 8];
           [TU_masks, ~] = get_masks(VTURU_glm, 'TU', clusterFWEcorrect, extent, Num);
           TU_mask1 = TU_masks{TU_rois(1)};
           TU_mask2 = TU_masks{TU_rois(2)};

           TU_betas1 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask1);
           TU_betas2 = get_beta_series(EXPT, beta_series_glm, subj, event, TU_mask2);

           multi.pmod(1).name{1} = 'RU_betas';
           multi.pmod(1).param{1} = RU_betas';
           multi.pmod(1).poly{1} = 1; 


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));


        % 36 but minimalist -- just trial onset, RU, and feedback onset
        % looking for RU ...
        %
        case 56 % <--- cluster in RLPFC is 230 voxels yay; still n.s....
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'feedback_onset';
           multi.onsets{2} = data(subj).feedback_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));


        % 23 but w/o choice onset => b/c WTF trial_onset betas from 23 in R RLPFC don't correlate with RU...
        %
        case 57 % <-- meh; going back to 23 b/c 58 and 59 kinda sucked (RU esp.)
           idx = 0;

           block = data(subj).block(which_trials);
           trial = data(subj).trial(which_trials);

           onsets = data(subj).trial_onset(which_trials);
           for t = 1:numel(onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_block_', num2str(block(t)), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['trial_onset_', suffix];
               multi.onsets{idx} = [onsets(t)];
               multi.durations{idx} = [0];
           end

           onsets = data(subj).feedback_onset(which_trials);
           for t = 1:numel(onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_block_', num2str(block(t)), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['feedback_onset_', suffix];
               multi.onsets{idx} = [onsets(t)];
               multi.durations{idx} = [0];
           end


        % 36 but w/o choice_onset
        % see 56
        %
        case 58  % <-- similar to 36; RU bigger RLPFC but wrong peak => no univariate_decoder; more or less same results for TU => reverting to 36 & 47
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'feedback_onset';
           multi.onsets{2} = data(subj).feedback_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'trial_onset_L';
           multi.onsets{3} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{3} = zeros(size(multi.onsets{3}));

        % 47 but w/o choice_onset
        %
        case 59 % <-- more or less the same as  47
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials); 
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'feedback_onset';
           multi.onsets{2} = data(subj).feedback_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'trial_onset_L';
           multi.onsets{3} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{3} = zeros(size(multi.onsets{3}));

        % 36 with decTU and decRU => functional connectivity
        %
        case 60  % <-- nothing for decRU and decTU !!!!! w t f well kinda makes sense but not really => shoulda been getter than RU and TU
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           data_cur = data; % we overwrite it temporarily
           load univariate_decoder_both_for_glmodel60_nosign.mat;

           decRU = data(subj).act(which_trials,1);
           decTU = data(subj).act(which_trials,2);

           multi.pmod(1).name{5} = 'decRU';
           multi.pmod(1).param{5} = decRU';
           multi.pmod(1).poly{5} = 1; 

           multi.pmod(1).name{6} = 'decTU';
           multi.pmod(1).param{6} = decTU';
           multi.pmod(1).poly{6} = 1; 

           data = data_cur;


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % #SAM
        % 47 with decRU and decTU => functional connectivity => do BIC / BMS thing in DV ROI from 47
        %
        case 61 % <-- nothing for RU nor TU; BUT BMS favors over 47
           [~, RU, TU, ~, DV] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials); 
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'DV';
           multi.pmod(1).param{1} = DV';
           multi.pmod(1).poly{1} = 1;    


           data_cur = data; % we overwrite it temporarily
           load univariate_decoder_both_for_glmodel60_nosign.mat;

           decRU = data(subj).act(which_trials,1);
           decTU = data(subj).act(which_trials,2);

           multi.pmod(1).name{2} = 'decRU';
           multi.pmod(1).param{2} = decRU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'decTU';
           multi.pmod(1).param{3} = decTU';
           multi.pmod(1).poly{3} = 1; 

           data = data_cur;


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % 36 with RPE
        %
        case 62
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');
           [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, RPE] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials & ~data(subj).timeout);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.pmod(3).name{1} = 'RPE';
           multi.pmod(3).param{1} = RPE';
           multi.pmod(3).poly{1} = 1;    

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'feedback_onset_timeouts';
               multi.onsets{5} = data(subj).feedback_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end

        % 36 with PRPE and NRPE
        %
        case 63
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');
           [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, RPE] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           feedback_onset = data(subj).feedback_onset(which_trials & ~data(subj).timeout);

           multi.names{3} = 'feedback_onset_pos';
           multi.onsets{3} = feedback_onset(RPE >= 0);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.pmod(3).name{1} = 'PRPE';
           multi.pmod(3).param{1} = RPE(RPE >= 0)';
           multi.pmod(3).poly{1} = 1;    

           multi.names{4} = 'feedback_onset_neg';
           multi.onsets{4} = feedback_onset(RPE < 0);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           multi.pmod(4).name{1} = 'NRPE';
           multi.pmod(4).param{1} = RPE(RPE < 0)';
           multi.pmod(4).poly{1} = 1;    

           multi.names{5} = 'trial_onset_L';
           multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{5} = zeros(size(multi.onsets{5}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = data(subj).feedback_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
           end


        % 36 without V/TU b/c that's pointless
        %
        case 64 % <-- same as 36
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % 64 with 1/TU
        %
        case 65 % <-- nothing for 1/TU
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'invTU';
           multi.pmod(1).param{4} = 1./TU';
           multi.pmod(1).poly{4} = 1; 


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

        % same as 45 but w/o choice_onset 
        %
        case 66
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'TU';
           multi.pmod(1).param{2} = TU';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V';
           multi.pmod(1).poly{3} = 1; 

           multi.pmod(1).name{4} = 'VTU';
           multi.pmod(1).param{4} = VTU';
           multi.pmod(1).poly{4} = 1; 

           multi.names{2} = 'feedback_onset';
           multi.onsets{2} = data(subj).feedback_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'trial_onset_L';
           multi.onsets{3} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{4} = 'trial_onset_timeouts';
               multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
           end



        % copy of 45 but with RU only -- for VIFs followup 
        % same idea as 39
        %
        case 67
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % copy of 45 but with TU only -- for VIFs followup 
        % same idea as 40
        %
        case 68
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'TU';
           multi.pmod(1).param{1} = TU';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % copy of 45 but with TU only -- for VIFs followup 
        % same idea as 41
        %
        case 69
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'V';
           multi.pmod(1).param{1} = V';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end


        % copy of 45 but with TU only -- for VIFs followup 
        % same idea as 42
        %
        case 70
        
           [V, RU, TU, VTU] = get_latents(data, subj, which_trials & ~data(subj).timeout, 'abs');

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials & ~data(subj).timeout);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'VTU';
           multi.pmod(1).param{1} = VTU';
           multi.pmod(1).poly{1} = 1;    

           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials & ~data(subj).timeout);
           multi.durations{2} = zeros(size(multi.onsets{2}));

           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).trial_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           if sum(which_trials & data(subj).timeout) > 0
               multi.names{5} = 'trial_onset_timeouts';
               multi.onsets{5} = data(subj).trial_onset(which_trials & data(subj).timeout); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
           end



        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end



