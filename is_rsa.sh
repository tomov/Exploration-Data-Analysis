mkdir output

outfileprefix="output/is_rsa"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running is_rsa  >> jobs.txt
echo ---------------- >> jobs.txt

#function is_rsa(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept)

declare -a fn_calls=(
                    "is_rsa"
                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo is_rsa.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

