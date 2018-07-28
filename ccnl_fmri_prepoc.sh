# run single subject GLM for a range of models for all subjects
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go ccnl_fmri_preproc >> jobs.txt
echo ---------------- >> jobs.txt

# this is crucial -- we can't simulate the subjects that don't have the full data; we don't want ccnl_fmri_glm to error out (note that when we were doing them in parallel before, it didn't matter if one subject failed b/c all other jobs would still continue)
#
goodSubjects=( 1 )  # getGoodSubjects()
subj_arg="${goodSubjects[@]}" # stringify it

outfileprefix="output/ccnl_fmri_preproc"
echo File prefix = $outfileprefix

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out --wrap="matlab -nodisplay -nosplash -nojvm -r $'ccnl_fmri_preproc(exploration_expt(),  [$subj_arg]);exit'"`
# for local testing
#sbatch_output=`echo Submitted batch job 88725418`
echo $sbatch_output

# Append job id to jobs.txt
#
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo ccnl_fmri_preproc.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
