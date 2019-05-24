% conjunction analysis

clear all;

EXPT = exploration_expt;

p = 0.001;
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3
clusterFWEcorrect = true;
extent = [];


con(1).glmodel = 45;
con(1).contrast = 'RU';
con(1).direct = '-';

con(2).glmodel = 45;
con(2).contrast = 'TU';
con(2).direct = '+/-';

%{
con(3).glmodel = 41;
con(3).contrast = 'V';
con(3).direct = '-';

con(4).glmodel = 47;
con(4).contrast = 'DV';
con(4).direct = '-';
%}


for c = 1:length(con)

    [con(c).V, con(c).Y, con(c).C, con(c).CI, con(c).region, con(c).extent, con(c).stat, con(c).mni, con(c).cor, con(c).results_table] = ccnl_extract_clusters(EXPT, con(c).glmodel, con(c).contrast, p, con(c).direct, alpha, Dis, Num, clusterFWEcorrect, extent);
end


conj = zeros(size(con(1).Y));

for c = 1:length(con)
    conj = conj + (con(c).CI > 0);
end


V = con(1).V;
V.fname = 'masks/conj_45.nii'; % change immediately!

spm_write_vol(V, conj);

struc = fullfile(EXPT.modeldir, 'mean.nii');
bspmview(V.fname, struc);

