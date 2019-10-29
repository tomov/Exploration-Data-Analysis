mkdir output

outfileprefix="output/neurosynth_CV"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running neurosynth_CV  >> jobs.txt
echo ---------------- >> jobs.txt


#function neurosynth_CV(regressor, do_orth, standardize, mixed_effects, intercept, method, get_null, zscore_across_voxels, predict_abs, use_smooth, lateralized, parcel_idxs)

declare -a fn_calls=(
                     "neurosynth_CV(\'RU\', false, 0, true, true, {\'fitrlinear_CV_1\', \'fitrlinear_CV_2\'}, false, false, false, false, true)"
                     "neurosynth_CV(\'TU\', false, 0, true, true, {\'fitrlinear_CV_1\', \'fitrlinear_CV_2\'}, false, false, false, false, true)"
                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 20-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo neurosynth_CV.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

