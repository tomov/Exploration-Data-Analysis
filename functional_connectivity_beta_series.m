% check for functional connectivity between regions by correlating their beta series (i.e. deconvolved activations) 
% WARNING: they will obviously be significant b/c we selected them to be so; we can only look at relative, trial onset vs feedback onset
% TODO copy of functional_connectivity_activations; dedupe

data = load_data;
[~,~,goodRuns] = exploration_getSubjectsDirsAndRuns();
EXPT = exploration_expt;
nTRs = 242;

filename = 'functional_connectivity_beta_series.mat'

ts = [];
dfs = [];
ps = [];
ts_to = [];
dfs_to = [];
ps_to = [];
ts_fo = [];
dfs_fo = [];
ps_fo = [];

reg_i = {};
reg_j = {};

clear corr_r;
clear corr_r_to;
clear corr_r_fo;
clear corr_z;
clear corr_z_to;
clear corr_z_fo;
clear region;

region(1).name = 'DV';
region(1).mask = 'sphere_glm29_DV_-38_-8_62_r=10mm.nii';
region(1).glm = 29;

region(2).name = 'RU';
region(2).mask = 'masks/badre_rlpfc_36_56_-8_r=10.0mm.nii';
region(2).glm = 45;

region(3).name = 'TU';
region(3).mask = 'masks/badre_dlpfc_38_30_34_r=10.0mm.nii';
region(3).glm = 45;

region(4).name = 'V';
region(4).mask = 'masks/vmpfc_l_resized.nii';
region(4).glm = 69;

for s = 1:length(data)
    % figure out which trials to exclude from hybrid GLMs (timeouts & bad runs with too much motion)
    runs = find(goodRuns{s}); % only those runs were included in the GLMs
    data(s).bad_runs = ~ismember(data(s).run, runs); % ... those runs were NOT included in the GLMs

    % figure out which TRs toughly correspond to trial_onsets (with HRF offset & stuff)
    %
    % trial onset idx's from bad runs are NaNs
    % note that timeouts are included here b/c we do have them in the fMRI GLM (not the behavioral one though)
    [~,session] = find(data(s).run == runs); % automatically excludes bad runs
    data(s).trial_onset_act_idx = nan(length(data(s).trial_onset), 1);
    data(s).trial_onset_act_idx(~data(s).bad_runs) = get_activations_idx(EXPT, data(s).trial_onset(~data(s).bad_runs), session, nTRs);

    % same for feedback onset
    data(s).feedback_onset_act_idx = nan(length(data(s).feedback_onset), 1);
    data(s).feedback_onset_act_idx(~data(s).bad_runs) = get_activations_idx(EXPT, data(s).feedback_onset(~data(s).bad_runs), session, nTRs);
end


beta_series_glm = 23;

for i = 1:length(region)
    for j = 1:i-1
        for s = 1:length(data)
            %act_i = {rand(nTRs*8)}; % ccnl_get_activations(EXPT, region(i).glm, region(i).mask, s);
            %act_j = {rand(nTRs*8)}; %ccnl_get_activations(EXPT, region(j).glm, region(j).mask, s);
            %act_i = ccnl_get_activations(EXPT, region(i).glm, region(i).mask, s);
            %act_j = ccnl_get_activations(EXPT, region(j).glm, region(j).mask, s);
            act_i_to = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', region(i).mask);
            act_i_fo = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'feedback_onset', region(i).mask);
            act_j_to = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'trial_onset', region(j).mask);
            act_j_fo = ccnl_get_beta_series(EXPT, beta_series_glm, s, 'feedback_onset', region(j).mask);
            act_i_to = nanmean(act_i_to, 2);
            act_i_fo = nanmean(act_i_fo, 2);
            act_j_to = nanmean(act_j_to, 2);
            act_j_fo = nanmean(act_j_fo, 2);

            % of trial onsets
            [r_to,p_to] = corr(act_i_to, act_j_to);
            corr_r_to{i,j}(s) = r_to;

            % of feedback onsets
            [r_fo,p_fo] = corr(act_i_fo, act_j_fo);
            corr_r_fo{i,j}(s) = r_fo;
        end

        % this is meaningless
        %{
        corr_z{i,j} = atanh(corr_r{i,j});
        [h,p,ci,stat] = ttest(corr_z{i,j});
        ts = [ts; stat.tstat];
        ps = [ps; p];
        dfs = [dfs; stat.df];
        %}

        corr_z_to{i,j} = atanh(corr_r_to{i,j});
        corr_z_fo{i,j} = atanh(corr_r_fo{i,j});


        [h,p,ci,stat] = ttest(corr_z_to{i,j}, corr_z_fo{i,j});
        ts = [ts; stat.tstat];
        ps = [ps; p];
        dfs = [dfs; stat.df];

        [h,p,ci,stat] = ttest(corr_z_to{i,j});
        ts_to = [ts_to; stat.tstat];
        ps_to = [ps_to; p];
        dfs_to = [dfs_to; stat.df];

        [h,p,ci,stat] = ttest(corr_z_fo{i,j});
        ts_fo = [ts_fo; stat.tstat];
        ps_fo = [ps_fo; p];
        dfs_fo = [dfs_fo; stat.df];


        reg_i = [reg_i; {region(i).name}];
        reg_j = [reg_j; {region(j).name}];
    end
end

save(filename, '-v7.3');

disp('trial onset - feedback onset');
table(reg_i, reg_j, ts, dfs, ps)

disp('trial onset');
table(reg_i, reg_j, ts_to, dfs_to, ps_to)

disp('feedback onset');
table(reg_i, reg_j, ts_fo, dfs_fo, ps_fo)
