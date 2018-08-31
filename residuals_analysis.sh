mkdir output

outfileprefix="output/residuals_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running residuals_analysis  >> jobs.txt
echo ---------------- >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(10, \'RU\', \'RU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 10 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'RU\', \'RU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'TU\', \'TU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
