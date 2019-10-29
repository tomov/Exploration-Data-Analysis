mkdir output

outfileprefix="output/multilinear_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running multilinear_analysis for TU - trial, GLM 21  >> jobs.txt
echo ---------------- >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'fitlm\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh fitlm : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'fitrlinear_ridge\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh fitrlinear_ridge: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'fitrlinear_ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh fitrlinear_ridge_CV: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'fitrlinear_lasso\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh fitrlinear_lasso: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'fitrlinear_lasso_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh fitrlinear_lasso_CV: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'ridge\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh ridge : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'ridge_CV\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo multilinear_analysis.sh ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'lasso\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh lasso : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1;
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'lasso_CV\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo multilinear_analysis.sh lasso_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
