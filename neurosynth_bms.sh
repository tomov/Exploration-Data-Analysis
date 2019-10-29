mkdir output

outfileprefix="output/neurosynth_bms"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running neurosynth_bms  >> jobs.txt
echo ---------------- >> jobs.txt


#function neurosynth_bms(regressor, do_orth, standardize, mixed_effects, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth, lateralized, parcel_idxs)

declare -a fn_calls=(
                     "neurosynth_bms(\'RU\', false, 0, true, true, \'fitrlinear_ridge\', false, false, false, false, true)"
                     "neurosynth_bms(\'TU\', false, 0, true, true, \'fitrlinear_ridge\', false, false, false, false, true)"
                     )

# double check
#declare -a fn_calls=(
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV_CV\', false, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge_CV_CV\', false, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV_CV\', true, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge_CV_CV\', true, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge\', false, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge\', false, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV\', false, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge_CV\', false, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge\', true, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge\', true, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV\', true, false, false, false, true, [18 18+200])"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge_CV\', true, false, false, false, true, [11 115 122 145 152  11+200 115+200 122+200 145+200 152+200])"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(\'RU\', false, 0, false, false, \'ridge\', false, false, false, false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, false, \'ridge\', false, false, false, false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge\', false, false, false, false)"
#                     "neurosynth_bms(\'TU\', false, 0, false, false, \'ridge\', false, false, false, false)"
#                     "neurosynth_bms(\'TU\', false, 0, true, false, \'ridge\', false, false, false, false)"
#                     "neurosynth_bms(\'TU\', false, 0, true, true, \'ridge\', false, false, false, false)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(\'RU\', false, 0, false, false, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, false, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, false, false, \'fitlm\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, false, false, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, false, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, true, true, \'ridge_CV\', false)"
#                     "neurosynth_bms(\'RU\', false, 0, false, false, \'fitlm\', false)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(47, \'DV\', 47, \'DV\', true, 1, 2, true, false)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'badre\', 36, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, true, false)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(47, \'DV\', 47, \'DV\', true, 1, 2, true, false)"
#                     "neurosynth_bms(47, \'DV\', 47, \'DV\', false, 1, 2, true, false)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', true, 1, 2, false, false)"
#                     )



#declare -a fn_calls=(
#                     "neurosynth_bms(45, \'RU\', 45, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(45, \'RU\', 45, \'TU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(45, \'TU\', 45, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(45, \'TU\', 45, \'TU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'badre\', 45, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'badre\', 45, \'TU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'dlpfc\', 45, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(-1, \'dlpfc\', 45, \'TU\', true, 1, 2, true, false)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true, false)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true, false)"
#                     )

# repro from paper
#declare -a fn_calls=(
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(36, \'RU\', 36, \'TU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(36, \'TU\', 36, \'RU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(-1, \'badre\', 36, \'RU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(-1, \'badre\', 36, \'TU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(-1, \'dlpfc\', 36, \'RU\', true, 1, 2, false, false)"
#                     "neurosynth_bms(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, false, false)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(41, \'V\', 41, \'V\', false, 1, 2, true, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', false, 1, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', false, 1, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', false, 1, 2, true, true)"
#                     "neurosynth_bms(41, \'V\', 41, \'V\', false, 10, 2, true, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', false, 10, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', false, 10, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', false, 10, 2, true, true)"
#                     "neurosynth_bms(41, \'V\', 41, \'V\', false, 100, 2, true, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', false, 100, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', false, 100, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', false, 100, 2, true, true)"
#                     "neurosynth_bms(41, \'V\', 41, \'V\', false, 1000, 2, true, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', false, 1000, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', false, 1000, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', false, 1000, 2, true, true)"
#                     "neurosynth_bms(41, \'V\', 41, \'V\', false, 10000, 2, true, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', false, 10000, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', false, 10000, 2, true, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', false, 10000, 2, true, true)"
#                     )


#declare -a fn_calls=(
#                     "neurosynth_bms(41, \'V\', 41, \'V\', true, 1, 2, true)"
#                     "neurosynth_bms(41, \'V\', 36, \'RU\', true, 1, 2, true)"
#                     "neurosynth_bms(41, \'V\', 36, \'TU\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'RU\', 41, \'V\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'RU\', 36, \'TU\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'TU\', 41, \'V\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'RU\', true, 1, 2, true)"
#                     "neurosynth_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true)"
#                     )

#declare -a fn_calls=(
#                     "neurosynth_bms(11, \'RU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "neurosynth_bms(11, \'RU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     "neurosynth_bms(11, \'TU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "neurosynth_bms(11, \'TU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 20-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo neurosynth_bms.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

