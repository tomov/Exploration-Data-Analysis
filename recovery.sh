# run recovery for a bunch of subjects; must edit run_recovery.m first (that's where the action is)
#

mkdir output

outfileprefix="output/recovery"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running recovery  >> jobs.txt
echo ---------------- >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf_holy --mem 50000 -t 10-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'recovery;exit'"`
# for local testing
#sbatch_output=`echo Submitted batch job 88725418`
echo $sbatch_output

# Append job id to jobs.txt
#
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo recovery.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
