mkdir output

outfileprefix="output/cross_subj_lik"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subj_lik >> jobs.txt
echo ---------------- >> jobs.txt

#function cross_subj_lik(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)


declare -a fn_calls=(
                     "cross_subj_lik(29, \'DV\', 29, \'DV\', 2, false, 100, false)"
                     "cross_subj_lik(29, \'DV\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subj_lik(29, \'DV\', 45, \'TU\', 2, false, 100, false)"
                     "cross_subj_lik(29, \'DV\', 45, \'V\', 2, false, 100, false)"
                     "cross_subj_lik(-1, \'badre\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subj_lik(-1, \'badre\', 45, \'TU\', 2, false, 100, false)"
                     "cross_subj_lik(-1, \'dlpfc\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subj_lik(-1, \'dlpfc\', 45, \'TU\', 2, false, 100, false)"
                     )

# old
#declare -a fn_calls=(
#                     "cross_subj_lik(47, \'DV\', 47, \'DV\', 2, false, 100, false)"
#                     "cross_subj_lik(36, \'RU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subj_lik(36, \'RU\', 36, \'TU\', 2, false, 100, false)"
#                     "cross_subj_lik(36, \'TU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subj_lik(36, \'TU\', 36, \'TU\', 2, false, 100, false)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo cross_subj_lik.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

