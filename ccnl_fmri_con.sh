# run ccnl_fmri_con for a bunch of subjects; must edit run_ccnl_fmri_con.m first (that's where the action is)
#

mkdir output

outfileprefix="output/ccnl_fmri_con"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running ccnl_fmri_con  >> jobs.txt
echo ---------------- >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf_holy --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'run_ccnl_fmri_con;exit'"`
# for local testing
#sbatch_output=`echo Submitted batch job 88725418`
echo $sbatch_output

# Append job id to jobs.txt
#
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo ccnl_fmri_con.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
