% Pull beta series from a seed region
%
% adapted from https://github.com/hayleydorfman/agency-fmri/blob/master/getBetaSeries.m
%
function B = get_beta_series(EXPT,glmodel,subj, name, mask)

    %prefix = [event_prefix, '_run_', num2str(run)];
    
    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');

    % load betas 
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

    which = contains(SPM.xX.name, name); % betas for given event
    %which(which) = rsa.which_betas; % of those, only betas for given trials
    cdir = pwd;
    cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
    B = spm_data_read(SPM.Vbeta(which), find(mask));
    cd(cdir);
            
    B = nanmean(B,2);
    
    assert(size(B,1) > 0, 'no betas - likely wrong name');
