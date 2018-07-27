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
    
    [allSubjects, subjdirs, nRuns] = exploration_getSubjectsDirsAndRuns();
    %assert(isequal(allSubjects, metadata.allSubjects));
        
    % pick the trials that correspond to that subject & run
    % notice that we're not &-ing with data.which_rows -- we might want to
    % run the GLM even for subjects that are not good. In fact, we cannot
    % determine whether subjects are good or not before we have run the GLM
    % and inspected for stuff like rotation and such.
    %
    which_trials = data.subject == subj & data.run == run;
    assert(sum(which_trials) == 40);
    
    % ...never mind what the thing above said;
    % we only support good subjects here now
    %
    %assert(ismember(subj, getGoodSubjects()));
    
    
    
    % Parametric modulators
    %
    switch glmodel
        
        % Simple. Contrast conditions at trial onset 
        %
        case 1 % <------------- GOOD
            % condition @ trial onset 
            % 
            for cond = 1:4
                multi.names{cond} = conds{cond};
                multi.onsets{cond} = data.trial_onset(which_trials & data.cond == cond)';
                multi.durations{cond} = zeros(size(multi.onsets{cond}));
            end
            
            % nuisance @ choice onset
            % 
            multi.names{5} = 'choice_onset';
            multi.onsets{5} = data.choice_onset(which_trials)';
            multi.durations{5} = zeros(size(multi.onsets{5}));

            % nuisance @ feedback onset
            % 
            multi.names{6} = 'feedback_onset';
            multi.onsets{6} = data.feedback_onset(which_trials)';
            multi.durations{6} = zeros(size(multi.onsets{6}));
        

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end 
