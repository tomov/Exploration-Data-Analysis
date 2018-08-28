function [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx)

    % Computes behavioral RDMs.
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % EXAMPLE:
    %   [Behavioral, control] = ccnl_behavioral_rdms(exploration_expt(), 1);
    %   showRDMs(Behavioral(1).subj);
    %
    % INPUT:
    %   EXPT = experiment structure
    %   rsa_idx - which RSA to use 
    %
    % OUTPUTS:
    %   Behavioral - struct array, one element per model, with the fields:
    %      .name - model name (e.g. 'condition')
    %      .D - dimension of features
    %      .is_control - boolean whether this is a control model
    %      .distance_measure - name (e.g. 'cosine') or function handler to be used as a distance measure for the RDMs (passed to pdist)
    %      .subj - struct array with subject RDMs for given model; has following fields:
    %         .features - [nTrials x D] feature vector
    %         .RDM - [nTrials x nTrials] RDM based on features
    %   control - indices of models in Behavioral that are controls
    %
    
    for s = 1:length(EXPT.subject) % for each subject
        rsa = EXPT.create_rsa(rsa_idx, s, 1);
        n = length(rsa.model);

        % initialize model metadata for all models
        for i = 1:n % for each RSA model
            if s == 1
                Behavioral(i).D = size(rsa.model(i).features, 2);
                Behavioral(i).name = rsa.model(i).name;
                Behavioral(i).is_control = rsa.model(i).is_control;
                Behavioral(i).distance_measure = rsa.model(i).distance_measure;
            else
                assert(Behavioral(i).D == size(rsa.model(i).features, 2));
                assert(isequal(Behavioral(i).name, rsa.model(i).name));
                assert(Behavioral(i).is_control == rsa.model(i).is_control);
            end
            Behavioral(i).subj(s).features = zeros(0, Behavioral(i).D);
            Behavioral(i).subj(s).name = Behavioral(i).name;
        end

        % set data for each model for all runs
        for r = 1:length(EXPT.subject(s).functional) % for each run
            rsa = EXPT.create_rsa(rsa_idx, s, r);
            for i = 1:n % for each model
                Behavioral(i).subj(s).features = [Behavioral(i).subj(s).features; rsa.model(i).features];
            end
        end
    end

    % compute RDMs
    for i = 1:length(Behavioral)
        for s = 1:length(Behavioral(i).subj)
            Behavioral(i).subj(s).RDM = squareRDMs(pdist(Behavioral(i).subj(s).features, Behavioral(i).distance_measure));
        end
    end

    control = find([Behavioral.is_control]);
end

