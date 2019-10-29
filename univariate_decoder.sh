mkdir output

outfileprefix="output/univariate_decoder"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running univariate_decoder  >> jobs.txt
echo ---------------- >> jobs.txt

#function univariate_decoder(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign)

declare -a fn_calls=(
                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(36, \'RU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(36, \'TU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(-1, \'badre\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(-1, \'badre\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(-1, \'dlpfc\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true)"
                     "univariate_decoder(-1, \'dlpfc\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true)"
                     )

#declare -a fn_calls=(
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, true)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 0.0001, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 0.001, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 0.01, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 0.1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 10, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 100, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 1000, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 10000, 2, false, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'badre\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(47, \'DV\', 47, \'DV\', false, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )



#declare -a fn_calls=(
#                     "univariate_decoder(45, \'RU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(45, \'RU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(45, \'TU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(45, \'TU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'badre\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'badre\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'dlpfc\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'dlpfc\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

# repro from paper
#declare -a fn_calls=(
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(36, \'RU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'badre\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'badre\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'dlpfc\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(41, \'V\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(41, \'V\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(41, \'V\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(41, \'V\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(41, \'V\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', false, 10000, 2, true, false, 100, 1, true)"
#                     )


#declare -a fn_calls=(
#                     "univariate_decoder(41, \'V\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(41, \'V\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(41, \'V\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'RU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'RU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'TU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'TU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder(11, \'RU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder(11, \'RU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder(11, \'TU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder(11, \'TU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo univariate_decoder.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

