mkdir output

outfileprefix="output/cross_subj_lik_bic"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subj_lik_bic >> jobs.txt
echo ---------------- >> jobs.txt

#function cross_subj_lik_bic(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs, contra)

# conjunction, GLMs 36 and 64
declare -a fn_calls=(
                     "cross_subj_lik_bic(29, \'DV\', 29, 0, false, 100, false, false)"
                     "cross_subj_lik_bic(29, \'DV\', 29, 0, false, 100, false, true)"
                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 0-3:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo cross_subj_lik_bic.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

