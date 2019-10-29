mkdir output

outfileprefix="output/univariate_decoder_both_badre"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running univariate_decoder_both_badre  >> jobs.txt
echo ---------------- >> jobs.txt

#function univariate_decoder_both(glmodel, RU_roi_idx, TU_roi_idx, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign)
declare -a fn_calls=(
                     "univariate_decoder_both_badre(45, 1, 2, false, 1, 2, true, false, 100, 1, true, true)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo univariate_decoder_both_badre.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

