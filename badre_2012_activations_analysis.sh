mkdir output

outfileprefix="output/badre_2012_activations_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running badre_2012_activations_analysis  >> jobs.txt
echo ---------------- >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_activations_analysis(21, 0);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_activations_analysis.sh normalize=0 glm 21 non-normalize: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt


sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_activations_analysis(21, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo badre_2012_activations_analysis.sh glm 21 normalize=1 no correction for abs value: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_activations_analysis(21, 2);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_activations_analysis.sh glm 21 normalize=2: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
