mkdir output

outfileprefix="output/main_effect"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running main_effect >> jobs.txt
echo ---------------- >> jobs.txt

# main_effect(roi_glmodel, roi_contrast, glmodel, regressor, clusterFWEcorrect, extent, Num)

#declare -a fn_calls=(
#                     "main_effect(47, \'DV\', 52, \'RU_betas\', false, 100, 1)"
#                     "main_effect(47, \'DV\', 52, \'TU_betas1\', false, 100, 1)"
#                     "main_effect(47, \'DV\', 52, \'TU_betas2\', false, 100, 1)"
#                     "main_effect(47, \'DV\', 53, \'RU_betas\', false, 100, 1)"
#                     "main_effect(47, \'DV\', 53, \'TU_betas1\', false, 100, 1)"
#                     "main_effect(47, \'DV\', 53, \'TU_betas2\', false, 100, 1)"
#                     )


# controls
#declare -a fn_calls=(
#                     "main_effect(-1, \'dlpfc\', 36, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 36, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 39, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 40, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 45, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 45, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 46, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 46, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 56, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 62, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 62, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 64, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 64, \'TU\', false, 100, 1)"
#                     )

# main 
declare -a fn_calls=(
                     "main_effect(-1, \'badre\', 36, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 36, \'TU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 39, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 40, \'TU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 45, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 45, \'TU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 46, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 46, \'TU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 56, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 62, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 62, \'TU\', false, 100, 1)"
                     "main_effect(-1, \'badre\', 64, \'RU\', false, 100, 1)"
                     "main_effect(-1, \'dlpfc\', 64, \'TU\', false, 100, 1)"
                     )

# GLM 66
#declare -a fn_calls=(
#                     "main_effect(-1, \'badre\', 66, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 66, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 66, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 66, \'TU\', false, 100, 1)"
#                     )

# GLM 45 diff spheres
#declare -a fn_calls=(
#                     "main_effect(-1, \'badre\', 45, \'RU\', false, 100, 1, 15)"
#                     "main_effect(-1, \'dlpfc\', 45, \'TU\', false, 100, 1, 15)"
#                     "main_effect(-1, \'badre\', 45, \'RU\', false, 100, 1, 5)"
#                     "main_effect(-1, \'dlpfc\', 45, \'TU\', false, 100, 1, 5)"
#                     "main_effect(-1, \'badre\', 45, \'RU\', false, 100, 1, 1)"
#                     "main_effect(-1, \'dlpfc\', 45, \'TU\', false, 100, 1, 1)"
#                     )

# looking for effects for badre ROIs
#declare -a fn_calls=(
#                     "main_effect(-1, \'badre\', 36, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 36, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 39, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 40, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 45, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 45, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 46, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 46, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 56, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 62, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 62, \'TU\', false, 100, 1)"
#                     "main_effect(-1, \'badre\', 64, \'RU\', false, 100, 1)"
#                     "main_effect(-1, \'dlpfc\', 64, \'TU\', false, 100, 1)"
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
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 0-3:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo main_effect.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done


