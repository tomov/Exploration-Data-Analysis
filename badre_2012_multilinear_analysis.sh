mkdir output

outfileprefix="output/badre_2012_multilinear_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running badre_2012_multilinear_analysis   >> jobs.txt
echo ---------------- >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_multilinear_analysis(\'ridge\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_multilinear_analysis.sh ridge : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_multilinear_analysis(\'ridge_CV\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo badre_2012_multilinear_analysis.sh ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_multilinear_analysis(\'lasso\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_multilinear_analysis.sh lasso : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_multilinear_analysis(\'lasso_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_multilinear_analysis.sh lasso_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'badre_2012_multilinear_analysis(\'fitrlinear_ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo badre_2012_multilinear_analysis.sh fitrlinear_ridge_CV: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
