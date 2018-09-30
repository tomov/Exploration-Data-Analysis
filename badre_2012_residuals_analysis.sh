mkdir output

outfileprefix="output/badre_2012_residuals_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running badre_2012_residuals_analysis  >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_residuals_analysis(21, false);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo badre_2012_residuals_analysis.sh glm 21 non-normalize: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt


sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_residuals_analysis(21, true);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo badre_2012_residuals_analysis.sh glm 21 normalize: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
