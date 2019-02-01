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
    %   multi - a structure with the folloowing fields
    %        .names{i}
    %        .onsets{i}
    %        .durations{i}
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
        case 36
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









        % 
        % ================= John's stuff ===============================
        %


        % =========== G_tot, N_tot =======================


        % ---------------------- fixed effects ---------------------------


        % #john
        % AU, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1001
           load('fit_AU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % ACU, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1002
           load('fit_ACU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % OpAL, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1003
           load('fit_OpAL_25nstarts_fixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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



        % ---------------------- mixed effects ---------------------------


        % #john
        % AU, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1004
           load('fit_AU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % ACU, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1005
           load('fit_ACU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % OpAL, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1006
           load('fit_OpAL_25nstarts_mixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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




        % ---------------------- random effects ---------------------------


        % #john
        % AU, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1007
           load('fit_AU_25nstarts_random.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % ACU, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1008
           load('fit_ACU_25nstarts_random.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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


        % #john
        % OpAL, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1009
           load('fit_OpAL_25nstarts_random.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'N';
           multi.pmod(1).param{2} = N';
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



        % =========== G_tot + N_tot =======================


        % ---------------------- fixed effects ---------------------------


        % #john
        % AU, fixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1010
           load('fit_AU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % ACU, fixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1011
           load('fit_ACU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % OpAL, fixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1012
           load('fit_OpAL_25nstarts_fixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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



        % ---------------------- mixed effects ---------------------------


        % #john
        % AU, mixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1013
           load('fit_AU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % ACU, mixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1014
           load('fit_ACU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % OpAL, mixed effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1015
           load('fit_OpAL_25nstarts_mixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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




        % ---------------------- random effects ---------------------------


        % #john
        % AU, random effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1016
           load('fit_AU_25nstarts_random.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % ACU, random effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1017
           load('fit_ACU_25nstarts_random.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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


        % #john
        % OpAL, random effects
        % G_tot + N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1018
           load('fit_OpAL_25nstarts_random.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'both';
           multi.pmod(1).param{1} = latents.a * G' + latents.b * N';
           multi.pmod(1).poly{1} = 1;

           multi.pmod(1).name{2} = 'diff';
           multi.pmod(1).param{2} = latents.a * G' - latents.b * N';
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



        % =========== G_tot =======================


        % ---------------------- fixed effects ---------------------------


        % #john
        % AU, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1019
           load('fit_AU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % ACU, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1020
           load('fit_ACU_25nstarts_fixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % OpAL, fixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1021
           load('fit_OpAL_25nstarts_fixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x);
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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



        % ---------------------- mixed effects ---------------------------


        % #john
        % AU, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1022
           load('fit_AU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % ACU, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1023
           load('fit_ACU_25nstarts_mixed.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % OpAL, mixed effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1024
           load('fit_OpAL_25nstarts_mixed.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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




        % ---------------------- random effects ---------------------------


        % #john
        % AU, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1025
           load('fit_AU_25nstarts_random.mat', 'results'); % only load results!

           latents = AU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % ACU, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1026
           load('fit_ACU_25nstarts_random.mat', 'results'); % only load results!

           latents = ACU(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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


        % #john
        % OpAL, random effects
        % G_tot, N_tot @ trial_onset
        % left choice @ trial_onset
        % nuisance @ choice_onset and feedback_onset
        %
        case 1027
           load('fit_OpAL_25nstarts_random.mat', 'results'); % only load results!

           latents = OpAL(data(subj), results.x(subj,:));
           G = sum(latents.G(which_trials,:), 2);
           N = sum(latents.N(which_trials,:), 2);

           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them

           multi.pmod(1).name{1} = 'G';
           multi.pmod(1).param{1} = G';
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














        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end



