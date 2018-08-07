
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
    
    

    [allSubjects, subjdirs, logical_vals] = exploration_getSubjectsDirsAndRuns();
    
    
    
    
    new_vals{subj} = ~logical_vals{subj};
    new_vals{subj} = cumsum(new_vals{subj});
    
    
    fprintf("%d",new_vals{subj});
    
     
    %run = run + new_vals{subj}(run); %change formula here
        
  
    fprintf('run %d \n', run);
    
    
    


    load(fullfile('results_glme_fig3.mat'), 'results_V', 'results_VTU', 'results_VRU', 'results_VTURU');

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
    





    %assert(isequal(allSubjects, metadata(subj).allSubjects));

    % pick the trials that correspond to that subject & run
    % notice that we're not &-ing with data(subj).which_rows -- we might want to
    % run the GLM even for subjects that are not good. In fact, we cannot
    % determine whether subjects are good or not before we have run the GLM
    % and inspected for stuff like rotation and such.
    %
    which_trials = data(subj).run == run;  
    %which_trials = run;
    %assert(sum(which_trials) == 40);
    
    fprintf('%d', which_trials);
    
    

    % ...never mind what the thing above said;
    % we only support good subjects here now
    %
    [~,~,~,goodSubjs] = exploration_getSubjectsDirsAndRuns();
    assert(ismember(subj, goodSubjs));
    
    [results_V, results_VTU, results_VRU, results_VTURU ] = model_comparison(data(subj));





    % Parametric modulators
    %
    switch glmodel

        % Condition @ trial_onset 
        % nuisance @ choice and feedback onset 
        
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
          
        case 2
            
           latents = kalman_filter(data(subj));
           
        
   
           Ws = fixedEffects(results_VTURU);
           
          
           
           QL = latents.m(:,1);
           QR = latents.m(:,2);
           stdL = latents.s(:,1);
           stdR = latents.s(:,2);
           
           TU = sqrt(stdL + stdR);
           
           L = (Ws(1)*QL) + (Ws(2)*stdL) + ((QL./TU)*Ws(3));  
           R = (Ws(1)*QR) + (Ws(2)*stdL) + ((QR./TU)*Ws(3));
           
           
          for i=1:length(L)
              if L(i) >= R(i)
                  RU(i) = stdL(i) - stdR(i); 
                  V(i) = QL(i) - QR(i);
                  
              else 
                  RU(i) = stdR(i) - stdL(i);
                  V(i) = QR(i) - QL(i);
              end 
          end
          
           
   
           
           %getMax = max(latents.s(:,1),latents.s(:,2));
           %getMin = max(latents.s(:,1),latents.s(:,2));
           
           
           %RU = sqrt(latents.s(:,1)) - sqrt(latents.s(:,2));
           %V = latents.m(:,1) - latents.m(:,2);
           VTU = V./TU;


           multi.names{1} = 'trial_onset';
           multi.onsets{1} = data(subj).trial_onset(which_trials);
           multi.durations{1} = zeros(size(multi.onsets{1}));

           multi.orth{1} = 0; % do not orthogonalise them  

           multi.pmod(1).name{1} = 'RU';
           multi.pmod(1).param{1} = RU(which_trials)';
           multi.pmod(1).poly{1} = 1;    

           multi.pmod(1).name{2} = 'VTU';
           multi.pmod(1).param{2} = VTU(which_trials)';
           multi.pmod(1).poly{2} = 1; 

           multi.pmod(1).name{3} = 'V';
           multi.pmod(1).param{3} = V(which_trials)';
           multi.pmod(1).poly{3} = 1; 


           multi.names{2} = 'choice_onset';
           multi.onsets{2} = data(subj).choice_onset(which_trials);
           multi.durations{2} = zeros(size(multi.onsets{2}));


           multi.names{3} = 'feedback_onset';
           multi.onsets{3} = data(subj).feedback_onset(which_trials);
           multi.durations{3} = zeros(size(multi.onsets{3}));

           multi.names{4} = 'trial_onset_L';
           multi.onsets{4} = data(subj).choice_onset(which_trials & data(subj).choice == 1);
           multi.durations{4} = zeros(size(multi.onsets{4}));

           multi.names{5} = 'trial_onset_R';
           multi.onsets{5} = data(subj).choice_onset(which_trials & data(subj).choice == 2);
           multi.durations{5} = zeros(size(multi.onsets{5}));

            
            
        % RU, V, V/TU @ trial_onset 
        % nuisance @ choice_onset and feedback_onset 
        
          
               

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end
