mkdir output

outfileprefix="output/activations_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running activations_analysis  >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'activations_analysis(19, \'RU\', \'RU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo activations_analysis.sh 19 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'activations_analysis(19, \'TU\', \'TU\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo activations_analysis.sh 19 TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1
