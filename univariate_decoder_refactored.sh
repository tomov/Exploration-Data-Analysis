mkdir output

outfileprefix="output/univariate_decoder_refactored"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running univariate_decoder_refactored  >> jobs.txt
echo ---------------- >> jobs.txt

#function univariate_decoder_refactored(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, do_CV, get_null)

# GLM 69, badre, mixed, flip_sign (as it should be), lambda = 1, intercept on
declare -a fn_calls=(
                     "univariate_decoder_refactored(69, \'V\', 69, \'V\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
                     )

# GLM 45, badre, mixed, flip_sign (as it should be), lambda = 1, intercept on
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo univariate_decoder_refactored.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

