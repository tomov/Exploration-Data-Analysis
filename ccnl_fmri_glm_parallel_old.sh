# run single subject GLM for a range of models for all subjects
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go ccnl_fmri_glm >> jobs.txt
echo ---------------- >> jobs.txt

for model in {3..2}
do
        for i in {1..3}
        do
            outfileprefix="output/ccnl_fmri_glm_${model}_subj_${i}"
            echo File prefix = $outfileprefix

            # send the job to NCF
            #
            sbatch_output=`sbatch -p ncf --mem 50000 -t 1-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $'ccnl_fmri_glm(exploration_expt(), $model, [$i]);exit'"`
            # for local testing
            #sbatch_output=`echo Submitted batch job 88725418`
            echo $sbatch_output

            # Append job id to jobs.txt
            #
            sbatch_output_split=($sbatch_output)
            job_id=${sbatch_output_split[3]}
            echo ccnl_fmri_glm.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
        done
done
