function rsa = exploration_create_rsa(rsa_idx, subj)

    % Create rsa structure, helper function for creating EXPT in
    % exploration_expt.m
    %
    % USAGE: rsa = exploration_create_rsa(model,subj)
    %
    % INPUTS:
    %   rsa_idx - positive integer indicating which RSA we're doing
    %   subj - integer specifying which subject is being analyzed
    %
    % OUTPUTS:
    %   rsa - a structure with the following fields:
    %     .glmodel - which GLM to use to get the trial-by-trial betas; make sure to have a unique regressor for each trial, e.g. 'trial_onset_1', 'trial_onset_2', etc.
    %     .event - which within-trial event to use for neural activity; used to pick the right betas (needs to be substring of the regressor name), e.g. 'trial_onset'
    %     .mask - path to .nii file, or 3D binary vector of voxels to consider
    %     .radius - searchlight radius in voxels
    %     .which_betas - logical mask for which betas (trials) front the GLM to include (e.g. not timeouts)
    %     .model - struct array describing the models used for behavioral RDMs (see Kriegeskorte et al. 2008) with the fields:
    %         .name - model name
    %         .features - [nTrials x D] feature vector
    %         .distance_measure - name (e.g. 'cosine') or function handler to be used as a distance measure for the RDMs (passed to pdist, see MATLAB documentation)
    %         .is_control - whether this is a control model (e.g. time)
    %
    % Momchil Tomov, Sep 2018


    fprintf('rsa %d, subj %d\n', rsa_idx, subj);

    data = load_data;
    conds = {'RS', 'SR', 'RR', 'SS'};
    

    [allSubjects, subjdirs, goodRuns, goodSubjs] = exploration_getSubjectsDirsAndRuns();
    
  
    % skip bad runs and timeouts
    runs = find(goodRuns{subj});
    bad_run = ~ismember(data(subj).run, runs);
    exclude = bad_run | data(subj).timeout;
    which_trials = ~exclude;
   
    fprintf('which_trials = %s\n', sprintf('%d', which_trials));

    % RSAs
    %
    switch rsa_idx

        % basic 
        %
        case 1
            rsa.event = 'trial_onset';
            rsa.glmodel = 23;
            rsa.radius = 10 / 1.5;
            rsa.mask = 'masks/mask.nii';
            rsa.which_betas = ~data(subj).timeout(~bad_run);

            %rsa.regressors = {'trial_onset_subj_1_run_1', ... etc.. };
            %rsa.radius = 2.6666;

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
            rsa.glmodel = 23;
            rsa.radius = 10 / 1.5;
            rsa.mask = 'masks/mask.nii';
            rsa.which_betas = ~data(subj).timeout(~bad_run);

            rsa.model(1).name = 'Qs';
            rsa.model(1).features = [Q1, Q2];
            rsa.model(1).features = rsa.model(1).features + rand(size(rsa.model(1).features)) * 0.0001; % no 0's
            rsa.model(1).distance_measure = 'cosine';
            rsa.model(1).is_control = false;

            rsa.model(2).name = 'sigmas';
            rsa.model(2).features = [std1, std2];
            rsa.model(2).features = rsa.model(2).features + rand(size(rsa.model(2).features)) * 0.0001; % no 0's
            rsa.model(2).distance_measure = 'cosine';
            rsa.model(2).is_control = false;

            % controls

            rsa.model(3).name = 'time';
            rsa.model(3).features = data(subj).trial_onset(which_trials);
            rsa.model(3).distance_measure = 'euclidean';
            rsa.model(3).is_control = true;

            rsa.model(4).name = 'run';
            rsa.model(4).features = data(subj).run(which_trials);
            rsa.model(4).distance_measure = @(c1, c2) c1 ~= c2;
            rsa.model(4).is_control = true;

        otherwise
            assert(false, 'invalid rsa_idx -- should be one of the above');

    end % end of switch statement

end
