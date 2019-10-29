mkdir output

outfileprefix="output/simulate_more"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running simulate_more  >> jobs.txt
echo ---------------- >> jobs.txt

#function simulate_more(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept)

declare -a fn_calls=(
                    "simulate_more(16, 10)"
                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo simulate_more.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

