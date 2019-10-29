mkdir output

outfileprefix="output/cross_subject"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subject >> jobs.txt
echo ---------------- >> jobs.txt

#function cross_subject(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)

declare -a fn_calls=(
                     "cross_subject(29, \'DV\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subject(29, \'DV\', 45, \'TU\', 2, false, 100, false)"
                     )

#declare -a fn_calls=(
#                     "cross_subject(-1, \'badre\', 45, \'RU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'dlpfc\', 45, \'TU\', 2, false, 100, false)"
#                     )

# conjunction, GLM 45, mean and peak 
#declare -a fn_calls=(
#                     "cross_subject(-1, \'conj_3\', 45, \'RU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 45, \'TU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 45, \'V\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 45, \'RU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_3\', 45, \'TU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_3\', 45, \'V\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_45\', 45, \'RU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_45\', 45, \'TU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_45\', 45, \'V\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_45\', 45, \'RU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_45\', 45, \'TU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_45\', 45, \'V\', 2, false, 100, false, true)"
#                     )

# conjunction, separate GLMs, with mean and peak beta
#declare -a fn_calls=(
#                     "cross_subject(-1, \'conj_3\', 39, \'RU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 40, \'TU\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 41, \'V\', 2, false, 100, false)"
#                     "cross_subject(-1, \'conj_3\', 39, \'RU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_3\', 40, \'TU\', 2, false, 100, false, true)"
#                     "cross_subject(-1, \'conj_3\', 41, \'V\', 2, false, 100, false, true)"
#                     )

# separate GLMs
#declare -a fn_calls=(
#                     "cross_subject(39, \'RU\', 39, \'RU\', 2, false, 100, false)"
#                     "cross_subject(40, \'TU\', 40, \'TU\', 2, false, 100, false)"
#                     "cross_subject(41, \'V\', 41, \'V\', 2, false, 100, false)"
#                     )

# functional connectivity
#declare -a fn_calls=(
#                     "cross_subject(47, \'DV\', 52, \'RU_betas\', 2, false, 100, false)"
#                     "cross_subject(47, \'DV\', 52, \'TU_betas1\', 2, false, 100, false)"
#                     "cross_subject(47, \'DV\', 52, \'TU_betas2\', 2, false, 100, false)"
#                     "cross_subject(47, \'DV\', 53, \'RU_betas\', 2, false, 100, false)"
#                     "cross_subject(47, \'DV\', 53, \'TU_betas1\', 2, false, 100, false)"
#                     "cross_subject(47, \'DV\', 53, \'TU_betas2\', 2, false, 100, false)"
#                     )

#declare -a fn_calls=(
#                     "cross_subject(36, \'RU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subject(36, \'RU\', 36, \'TU\', 2, false, 100, false)"
#                     "cross_subject(36, \'TU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subject(36, \'TU\', 36, \'TU\', 2, false, 100, false)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 0-2:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo cross_subject.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done


