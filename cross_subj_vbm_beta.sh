mkdir output

outfileprefix="output/cross_subj_vbm_beta"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subj_vbm_beta >> jobs.txt
echo ---------------- >> jobs.txt

declare -a fn_calls=(
                     "cross_subj_vbm_beta(41, \'V\', \'V\', 2, false, 100, 1)"
                     "cross_subj_vbm_beta(36, \'RU\', \'RU\', 2, false, 100, 1)"
                     "cross_subj_vbm_beta(36, \'TU\', \'TU\', 2, false, 100, 1)"
                     )
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm(21, \'TU\', \'TU - trial\');exit'"`


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo run $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done


##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'TU - trial\', 2);exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, TU, TU - trial, standardize=2: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
#
#
#
##
##
##echo -- RU -- >> jobs.txt
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'RU\', \'dlpfc\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, RU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'RU\', \'badre\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, RU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'RU\', \'tommy\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, RU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'RU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subj_vbm_beta.sh GLM 21, RU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'RU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subj_vbm_beta.sh GLM 21, RU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
##
##
##
##
##echo -- TU -- >> jobs.txt
##
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'dlpfc\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, TU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'badre\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, TU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'tommy\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_vbm_beta.sh GLM 21, TU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'RU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subj_vbm_beta.sh GLM 21, TU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_vbm_beta(21, \'TU\', \'TU - trial\');exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subj_vbm_beta.sh GLM 21, TU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#
