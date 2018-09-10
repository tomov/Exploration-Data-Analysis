
experiment="Exploration" # must match directory name

# list of subject ID's as registered on CBS, e.g.  subjects=('180725_UEP_001' '189725_UEP_002')
subjects=('180725_UEP_001' '180727_UEP_002' '180727_UEP_003'  '180730_UEP_004' '180730_UEP_005' '180801_UEP_006' '180802_UEP_007' '180803_UEP_008'  '180803_UEP_009' '180804_UEP010' '180804_UEP_011' '180804_UEP_012' '180804_UEP_013'  '180804_UEP_014' '180804_UEP_015' '180805_UEP_016' '180805_UEP_017' '180805_UEP_018_2','180805_UEP_019' '180805_UEP_020' '180805_UEP_021' '180806_UEP_022' '180806_UEP_023'  '180807_UEP_024' '180807_UEP_025' '180807_UEP_026' '180808_UEP_027' '180808_UEP_028'  '180808_UEP_029' '180809_UEP_030' '180809_UEP_031')

for subj in ${subjects[*]}; do
    cd /Volumes/MomchilfMRI/${experiment}/subjects/
    mkdir /Volumes/MomchilfMRI/${experiment}/subjects/${subj}/

    /Users/momchil/Dropbox/Research/exploration/code/ArcGet.py -a cbscentral -s ${subj} -r MEMPRAGE\ RMS,Minn_HCP_1.5mm_Task
done
