mkdir output

outfileprefix="output/cross_subject"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subject >> jobs.txt
echo ---------------- >> jobs.txt


sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'TU - trial\', 2, false, 100);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo cross_subject.sh GLM 21, TU, TU - trial, standardize=2, corr=0, extent=100: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'TU - trial\', 2);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, TU - trial, standardize=2: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1



#
#
#echo -- RU -- >> jobs.txt
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'RU\', \'dlpfc\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, RU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'RU\', \'badre\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, RU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'RU\', \'tommy\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, RU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'RU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, RU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'RU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, RU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
#
#echo -- TU -- >> jobs.txt
#
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'dlpfc\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'badre\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'tommy\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subject(21, \'TU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subject.sh GLM 21, TU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1



