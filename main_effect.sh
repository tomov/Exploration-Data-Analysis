mkdir output

outfileprefix="output/main_effect"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running main_effect >> jobs.txt
echo ---------------- >> jobs.txt

echo -- RU -- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'badre\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, RU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'tommy\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, RU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'RU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, RU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'TU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, RU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1




echo -- TU -- >> jobs.txt


sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'badre\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, TU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'tommy\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, TU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'RU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, TU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'TU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo main_effect.sh GLM 21, TU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1



