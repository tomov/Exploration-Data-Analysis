mkdir output

outfileprefix="output/cross_subj_perf_bic"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subj_perf_bic >> jobs.txt
echo ---------------- >> jobs.txt

#function cross_subj_perf_bic(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs, contra)

# conjunction, GLMs 36 and 64
declare -a fn_calls=(
                     "cross_subj_perf_bic(29, \'DV\', 29, 0, false, 100, false, false)"
                     "cross_subj_perf_bic(29, \'DV\', 29, 0, false, 100, false, true)"
                     )

#declare -a fn_calls=(
#                     "cross_subj_perf_bic(-1, \'conj_3\', 45, 0, false, 100, false)"
#                     "cross_subj_perf_bic(-1, \'conj_45\', 45, 0, false, 100, false)"
#                     )



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 0-3:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo cross_subj_perf_bic.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

