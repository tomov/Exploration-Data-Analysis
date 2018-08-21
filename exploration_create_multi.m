
function multi = exploration_create_multi(glmodel, subj, run, save_output)

    % Create multi structure, helper function for creating EXPT in
    % imageryExpt.m
    %
    % USAGE: multi = imagery_create_multi(model,subj,run)
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
    % Cody Kommers, July 2016

    if nargin < 4 || isempty(save_output)
        save_output = false;
    end

    fprintf('glm %d, subj %d, run %d\n', glmodel, subj, run);

    data = load_data;
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

    [V, RU, ~, VTU] = get_latents(data, subj, which_trials, 'left');



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
        case 11
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








        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end





% how = how to pick option 1 so as to disambiguate option 1 vs. option 2 and determine the sign of RU, V, etc.
%     how = left -> 1 = left, 2 = right
%     how = max -> 1 = max, 2 = min
%     how = chosen -> 1 = chosen, 2 = unchosen
%     how = abs -> return abs value; option varies depending on quantity
%
function [V, RU, TU, VTU, DV, DQ1, DQ2, Q1, Q2, std1, std2, DQL, DQR, QL, QR, stdL, stdR, w] = get_latents(data, subj, which_trials, how)
    latents = kalman_filter(data(subj));
   
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
    tbl = data2table(data);
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
