function [dec, dec_null] = decode_regressor_CV(EXPT, glmodel, regressor, mask, lambda, subjects, do_CV, null_iters, nTRs)

    % copied & modified from ccnl_decode_regressor; might update in ccnl-fmri later 

    % Invert the GLM (using ridge regression) to decode a regressor from the neural data.
    %
    % The formula is:
    %    dec = (y - sum x_i * b_i) * b_reg / (b_reg^2 + lambda)
    % where:
    %    dec = decoded regressor
    %    y = neural data
    %    x_i = regressor i
    %    b_i = beta coefficient for regressor i
    %    the sum is over all regressors i that are not our regressor
    %    b_reg = is the beta coefficient for our regressor
    %    lambda = ridge regularization parameter
    %
    % USAGE:
    %   dec = ccnl_decode_regressor(EXPT, glmodel, regressor, mask, lambda, [subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   glmodel - model number
    %   regressor - regressor name
    %   mask - a mask image name (e.g., 'mask.nii') in MNI or native space,
    %          a list of voxel indices in native space,
    %          a binary vector/mask in native space,
    %          or a list of voxels in native coordinates as a [N x 3] matrix
    %          (use mni2cor first if coordinates are in MNI space)
    %   lambda (optional) - regularization constant (defaults to 0, i.e. maximum likelihood estimate)
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   dec - {nSubjects} struct array with decoded regressor for each subject:
    %      dec{s} - [nTRs x nVoxels] decoded regressor for each TR for each voxel in the mask, for subject s
    %
    % EXAMPLE:
    %   dec = ccnl_decode_regressor(exploration_expt(), 21, 'RU', 'rlpfc.nii', 1)
 

    if ~exist('lambda', 'var')
        lambda = 0;
    end

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    % convert logicals to indices
    if strcmp(mask_format, 'mask') || islogical(mask)
        mask = find(mask);
    end

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

        % extract betas B and (filtered) design matrix X
        cdir = pwd;
        cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
        if strcmp(mask_format, 'cor')
            B = spm_data_read(SPM.Vbeta, 'xyz', mask');
        else
            B = spm_data_read(SPM.Vbeta, mask);
        end
        cd(cdir);
        X = SPM.xX.xKXs.X;

        % optional cross-validation:
        % set betas for each run = avg of betas from other runs
        %
        if do_CV
            B_CV = B;
            for i = 1:size(SPM.xX.name)
                name = SPM.xX.name{i};
                assert(startsWith(name, 'Sn ('));

                k = strfind(name, ' ');
                suffix = name(k(1)-1:end); % ') blabla*bf(1)'

                % find all versions of this regressor from the other runs 
                idx = find(endsWith(SPM.xX.name, suffix));
                idx(idx == i) = []; % exclude regressor from current run -- this is crucial, otherwise we're not CVing

                B_CV(i,:) = mean(B(idx,:), 1);
            end
            assert(immse(mean(B_CV,1), mean(B,1)) < 1e-9);
        end

        % separate our regressor from the rest
        names = SPM.xX.name'; % regressor names
        which_reg = contains(names, regressor);

        % separate X's and betas into matrices that do or don't have our regressor
        B_noreg = B(~which_reg, :);
        B_reg = B(which_reg, :);
        B_reg = repelem(B_reg, size(X, 1) / size(B_reg, 1), 1); % we need one for each TR b/c we're doing element-wise divison by b_RU
        X_noreg = X(:, ~which_reg);
        X_reg = X(:, which_reg);

        % extract activations, then whiten & apply filter (see spm_spm.m)
        act = ccnl_get_activations(EXPT, glmodel, mask, subj);
        act = spm_filter(SPM.xX.K,SPM.xX.W*act);

        % decode regressor
        dec{s} = (act - X_noreg * B_noreg) .* B_reg ./ (B_reg.^2 + lambda);
        dec{s} = mean(dec{s}, 2); % average across voxels

        % optionally generate null distribution
        %
        dec_null{s} = [];
        assert(mod(size(act,1), nTRs) == 0);
        nruns = size(X,1) / nTRs;
        for r = 1:nruns
            idx = (r - 1) * nTRs + 1 : r * nTRs;
            act(idx,:) = act(randperm(idx,:))

            dec_null{s}(:,i) = (act - X_noreg * B_noreg) .* B_reg ./ (B_reg.^2 + lambda);
            dec_null{s}(:,i) = mean(dec_null{s}(:,i), 2);
        end
    end
