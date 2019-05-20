% sanity checking univariate decoder assumptions
% does decoded RU correlate with RU regressor, even before subsampling/deconvolution?
%
% find(dec_ps < 0.05)   => 12 subjects
%  1     3     6    10    14    16    17    20    22    28    30    31

EXPT = exploration_expt();
glmodel = 36;
regressor = 'RU';
mask = 'sphere_glm36_RU_34_48_-8_r=1mm.nii';
lambda = 0;


% load mask
[mask_format, mask, Vmask] = get_mask_format_helper(mask);

% convert logicals to indices
if strcmp(mask_format, 'mask') || islogical(mask)
    mask = find(mask);
end
if ~exist('subjects', 'var')
    subjects = 1:length(EXPT.subject);
end


for s = 1:length(subjects)


    subj = subjects(s);
    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

    % extract betas B and design matrix X
    cdir = pwd;
    cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
    if strcmp(mask_format, 'cor')
        B = spm_data_read(SPM.Vbeta, 'xyz', mask');
    else
        B = spm_data_read(SPM.Vbeta, mask);
    end
    cd(cdir);

    %X = SPM.xX.xKXs.X; %  use high-pass filtered & whitened design matrix
    X = SPM.xX.X; % use unfilterd, unwhitened design matrix

    % extract activations 
    act = ccnl_get_activations(EXPT, glmodel, mask, subj);
    %act = spm_filter(SPM.xX.K,SPM.xX.W*act); % whiten & high-pass filter (see spm_spm.m)

    % override B
    %B_old = B;
    %B = glmfit(X,act, 'normal','constant', 'off');

    % separate our regressor from the rest
    names = SPM.xX.name'; % regressor names
    which_reg = contains(names, regressor);

    assert(sum(which_reg) == length(EXPT.subject(subj).functional), 'Number of regressors that match is different from number of runs -- maybe add x and ^ prefix and suffix?');

    % separate X's and betas into matrices that do or don't have our regressor
    B_noreg = B(~which_reg, :);
    B_reg = B(which_reg, :);


    % are betas positive for that subject? if not, use different subject
    B_reg
    mean(B_reg)
    [h,p,ci,stat] = ttest(B_reg);
    fprintf('B_reg mean = %f, tstat = %f, p = %f\n', mean(B_reg), stat.tstat, p);


    % separate X's and betas into matrices that do or don't have our regressor
    B_reg = repelem(B_reg, size(X, 1) / size(B_reg, 1), 1); % we need one for each TR b/c we're doing element-wise divison by b_RU
    X_noreg = X(:, ~which_reg);
    X_reg = X(:, which_reg);



    % is RU (from X, i.e. convolved with hrf) correlated w/ activation
    RU = sum(X_reg, 2);
    [r,p] = corr(RU, act);

    act_rs(s) = r;
    act_ps(s) = p;

    fprintf('corr RU act: r = %f, p = %f\n', r, p);


    nTRs = 242;
    assert(mod(size(act,1), nTRs) == 0);
    assert(size(act,2) == 1, 'only 1 voxel for now; otherwise, do mean');
    nruns = size(act,1) / nTRs;

    %{
    rs = []; ps = [];
    for i = 1:nruns
        st = (i-1)*nTRs + 1;
        en = i * nTRs;
        x = X_reg(st:en, i);
        [r, p] = corr(x, act(st:en,:));
        rs = [rs r];
        ps = [ps p];
    end

    rs
    ps
    %}


    % decode regressor
    dec = (act - X_noreg * B_noreg) .* B_reg ./ (B_reg.^2 + lambda);


    % is decoded regressor correlated with RU?
    [r, p] = corr(RU, dec);

    fprintf('corr dec & RU : r = %f, p = %f\n', r, p);

    RUs{s} = RU;
    decs{s} = dec;

    dec_rs(s) = r;
    dec_ps(s) = p;


end


%{
% repeat w/ fixed effects design matrix => nothing...
%
regressor = 'trial_onset*bf';

sn1 = find(contains(names, 'Sn(1)'));
for i = 1:length(sn1)
    suffix = names{sn1(i)}(5:end);
    which = contains(names, suffix);
    Xf(:,i) = sum(X(:,which), 2);
end


y = act;

bhat = glmfit(Xf,y, 'normal','constant', 'off');

reg = find(contains(names(sn1), regressor));

[r,p] = corr(dec, X(:,reg))
%}
