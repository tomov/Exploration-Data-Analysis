mkdir output

outfileprefix="output/neural_ts_to_behavioral_weights"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running neural_ts_to_behavioral_weights  >> jobs.txt
echo ---------------- >> jobs.txt

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'RU\', \'RU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'RU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'TU\', \'TU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'TU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'RU\', \'RU - TU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 RU - TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(19, \'TU\', \'TU - RU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 19 TU - RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
#
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'RU\', \'RU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 21 RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'RU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 21 RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'TU\', \'TU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 21 TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'TU\', \'TU - trial\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo neural_ts_to_behavioral_weights.sh 21 TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'RU\', \'RU - TU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 21 RU - TU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'neural_ts_to_behavioral_weights(21, \'TU\', \'TU - RU\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo neural_ts_to_behavioral_weights.sh 21 TU - RU: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

