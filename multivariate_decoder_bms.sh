mkdir output

outfileprefix="output/multivariate_decoder_bms"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running multivariate_decoder_bms  >> jobs.txt
echo ---------------- >> jobs.txt

#function multivariate_decoder_bms(roi_glmodel, roi_contrast, regressor, do_orth, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth)

declare -a fn_calls=(
                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
                     "multivariate_decoder_bms(36, \'TU\', \'TU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'TU\', \'TU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, false, false, 100, 1, false, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, false, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'TU\', \'TU\', false, 0, false, false, 100, 1, false, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'TU\', \'TU\', false, 0, true, false, 100, 1, false, \'ridge_CV\', false, false, false, false)"
#                     "multivariate_decoder_bms(36, \'TU\', \'TU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false, false, false, false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, false, false, 100, 1, false, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, false, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(36, \'RU\', \'RU\', false, 0, false, false, 100, 1, false, \'fitlm\', false)"
#                     "multivariate_decoder_bms(56, \'RU\', \'RU\', false, 0, false, false, 100, 1, false, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(56, \'RU\', \'RU\', false, 0, true, false, 100, 1, false, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(56, \'RU\', \'RU\', false, 0, true, false, 100, 1, true, \'ridge_CV\', false)"
#                     "multivariate_decoder_bms(56, \'RU\', \'RU\', false, 0, false, false, 100, 1, false, \'fitlm\', false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'badre\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(47, \'DV\', 47, \'DV\', false, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )



#declare -a fn_calls=(
#                     "multivariate_decoder_bms(45, \'RU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(45, \'RU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(45, \'TU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(45, \'TU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'badre\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'badre\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'dlpfc\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'dlpfc\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

# repro from paper
#declare -a fn_calls=(
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'badre\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'badre\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'dlpfc\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "multivariate_decoder_bms(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', false, 10, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', false, 10, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', false, 100, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', false, 100, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', false, 10000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', false, 10000, 2, true, false, 100, 1, true)"
#                     )


#declare -a fn_calls=(
#                     "multivariate_decoder_bms(41, \'V\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(41, \'V\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(41, \'V\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'RU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'RU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'TU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "multivariate_decoder_bms(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     )

#declare -a fn_calls=(
#                     "multivariate_decoder_bms(11, \'RU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "multivariate_decoder_bms(11, \'RU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     "multivariate_decoder_bms(11, \'TU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "multivariate_decoder_bms(11, \'TU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo multivariate_decoder_bms.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

