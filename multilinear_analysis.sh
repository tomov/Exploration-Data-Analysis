mkdir output

outfileprefix="output/multilinear_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running multilinear_analysis   >> jobs.txt
echo ---------------- >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'TU\', \'TU - trial\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh GLM 21, TU - trial, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(\'badre\', \'RU\', \'RU\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh Badre ROI, RU, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(\'tommy\', \'TU\', \'TU\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh Tommy ROIs, TU, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'RU\', \'RU - TU\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh GLM 21, RU - TU, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt


sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'RU\', \'RU\', \'ridge_CV\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo multilinear_analysis.sh GLM 21, RU, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'RU\', \'RU - trial\', \'ridge_CV\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo multilinear_analysis.sh GLM 21, RU - trial, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

## controls


#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(21, \'RU\', \'TU - trial\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh GLM 21, RU control, TU - trial, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'multilinear_analysis(\'badre\', \'TU\', \'TU\', \'ridge_CV\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo multilinear_analysis.sh Badre ROI, TU control, ridge_CV : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
