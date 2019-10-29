mkdir output

outfileprefix="output/residuals_analysis"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running residuals_analysis w/ 1 peak voxel, with adjustment for abs RU, max beta voxel  >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'RU\', \'RU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'RU\', \'RU - trial\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'TU\', \'TU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'TU\', \'TU - trial\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1



sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'RU\', \'RU - TU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 RU - TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1


sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(19, \'TU\', \'TU - RU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 19 TU - RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1





sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'RU\', \'RU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'RU\', \'RU - trial\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'TU\', \'TU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'TU\', \'TU - trial\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'RU\', \'RU - TU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 RU - TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1


sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'residuals_analysis(21, \'TU\', \'TU - RU\', \'sphere\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo residuals_analysis.sh 21 TU - RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

