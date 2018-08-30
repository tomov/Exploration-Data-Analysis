mkdir output

outfileprefix="output/residual_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running residual_analysis  >> jobs.txt
echo ---------------- >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residual_analysis(10, \'RU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residual_analysis.sh random: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residual_analysis(19, \'RU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residual_analysis.sh random hierarchical: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residual_analysis(19, \'TU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residual_analysis.sh fixed: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
