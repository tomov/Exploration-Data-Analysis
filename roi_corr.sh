mkdir output

outfileprefix="output/roi_corr"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running roi_corr  >> jobs.txt
echo ---------------- >> jobs.txt

#function roi_corr(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept)

declare -a fn_calls=(
                    "roi_corr"
                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo roi_corr.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

