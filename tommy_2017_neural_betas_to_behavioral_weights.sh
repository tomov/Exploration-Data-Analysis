mkdir output

outfileprefix="output/tommy_2017_neural_betas_to_behavioral_weights"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running tommy_2017_neural_betas_to_behavioral_weights GLM 21   >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'tommy_2017_neural_betas_to_behavioral_weights;exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo tommy_2017_neural_betas_to_behavioral_weights.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

