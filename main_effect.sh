mkdir output

outfileprefix="output/main_effect"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running main_effect >> jobs.txt
echo ---------------- >> jobs.txt

# main_effect(roi_glmodel, roi_contrast, glmodel, regressor, clusterFWEcorrect, extent, Num)

declare -a fn_calls=(
                     "main_effect(47, \'DV\', 52, \'RU_betas\', false, 100, 1)"
                     "main_effect(47, \'DV\', 52, \'TU_betas1\', false, 100, 1)"
                     "main_effect(47, \'DV\', 52, \'TU_betas2\', false, 100, 1)"
                     "main_effect(47, \'DV\', 53, \'RU_betas\', false, 100, 1)"
                     "main_effect(47, \'DV\', 53, \'TU_betas1\', false, 100, 1)"
                     "main_effect(47, \'DV\', 53, \'TU_betas2\', false, 100, 1)"
                     )



#declare -a fn_calls=(
#                     "main_effect(-1, \'badre\', 36, \'RU\', false, 100, 3)"
#                     "main_effect(-1, \'badre\', 36, \'TU\', false, 100, 3)"
#                     "main_effect(-1, \'dlpfc\', 36, \'RU\', false, 100, 3)"
#                     "main_effect(-1, \'dlpfc\', 36, \'TU\', false, 100, 3)"
#                     )

#declare -a fn_calls=(
#                     "main_effect(41, \'V\', 41, \'V\', false, 100, 3)"
#                     "main_effect(41, \'V\', 36, \'RU\', false, 100, 3)"
#                     "main_effect(41, \'V\', 36, \'TU\', false, 100, 3)"
#                     "main_effect(36, \'RU\', 41, \'V\', false, 100, 3)"
#                     "main_effect(36, \'RU\', 36, \'RU\', false, 100, 3)"
#                     "main_effect(36, \'RU\', 36, \'TU\', false, 100, 3)"
#                     "main_effect(36, \'TU\', 41, \'V\', false, 100, 3)"
#                     "main_effect(36, \'TU\', 36, \'RU\', false, 100, 3)"
#                     "main_effect(36, \'TU\', 36, \'TU\', false, 100, 3)"
#                     )

#declare -a fn_calls=(
#                     "main_effect(11, \'RU\', 35, \'RU\', false, 100)"
#                     "main_effect(11, \'RU\', 35, \'TU\', false, 100)"
#                     "main_effect(11, \'TU\', 35, \'RU\', false, 100)"
#                     "main_effect(11, \'TU\', 35, \'TU\', false, 100)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo main_effect.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done




#
# old stuff
#

#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'badre\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, RU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'tommy\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, RU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, RU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'RU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, RU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
#
#echo -- TU -- >> jobs.txt
#
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'badre\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, TU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'tommy\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, TU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, TU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'main_effect(21, \'TU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo main_effect.sh GLM 21, TU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
