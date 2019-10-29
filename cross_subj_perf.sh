mkdir output

outfileprefix="output/cross_subj_perf"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running cross_subj_perf >> jobs.txt
echo ---------------- >> jobs.txt

#function cross_subj_perf(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)
#declare -a fn_calls=(
#                     "cross_subj_perf(47, \'DV\', 47, \'DV\', 2, false, 100, false)"
#                     )

declare -a fn_calls=(
                     "cross_subj_perf(29, \'DV\', 29, \'DV\', 2, false, 100, false)"
                     "cross_subj_perf(29, \'DV\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subj_perf(29, \'DV\', 45, \'TU\', 2, false, 100, false)"
                     "cross_subj_perf(29, \'DV\', 45, \'V\', 2, false, 100, false)"
                     "cross_subj_perf(-1, \'badre\', 45, \'RU\', 2, false, 100, false)"
                     "cross_subj_perf(-1, \'dlpfc\', 45, \'TU\', 2, false, 100, false)"
                     )

# conjunction, GLM 36, mean and peak 
#declare -a fn_calls=(
#                     "cross_subj_perf(-1, \'conj_3\', 45, \'RU\', 2, false, 100, false)"
#                     "cross_subj_perf(-1, \'conj_3\', 45, \'TU\', 2, false, 100, false)"
#                     "cross_subj_perf(-1, \'conj_3\', 45, \'V\', 2, false, 100, false)"
#                     "cross_subj_perf(-1, \'conj_45\', 45, \'RU\', 2, false, 100, false)"
#                     "cross_subj_perf(-1, \'conj_45\', 45, \'TU\', 2, false, 100, false)"
#                     "cross_subj_perf(-1, \'conj_45\', 45, \'V\', 2, false, 100, false)"
#                     )


#function cross_subj_perf(roi_glmodel, roi_contrast, glmodel, regressor, standardize, clusterFWEcorrect, extent, odd_runs)
#declare -a fn_calls=(
#                     "cross_subj_perf(36, \'RU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subj_perf(36, \'RU\', 36, \'TU\', 2, false, 100, false)"
#                     "cross_subj_perf(36, \'TU\', 36, \'RU\', 2, false, 100, false)"
#                     "cross_subj_perf(36, \'TU\', 36, \'TU\', 2, false, 100, false)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo cross_subj_perf.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done




#sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'TU - trial\', 2, false, 100);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo cross_subj_perf.sh GLM 21, TU, TU - trial, standardize=2, corr=0, extent=100: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'TU - trial\', 2);exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, TU - trial, standardize=2: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
#
#
#
##
##
##echo -- RU -- >> jobs.txt
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'RU\', \'dlpfc\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, RU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'RU\', \'badre\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, RU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'RU\', \'tommy\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, RU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'RU\', \'RU - trial\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, RU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'RU\', \'TU - trial\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, RU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##
##
##
##echo -- TU -- >> jobs.txt
##
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'dlpfc\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, dlpfc: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'badre\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, badre: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'tommy\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, tommy: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'RU - trial\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, RU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
##
##sbatch_output=`sbatch -p ncf --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'cross_subj_perf(21, \'TU\', \'TU - trial\');exit'"`
##sbatch_output_split=($sbatch_output)
##job_id=${sbatch_output_split[3]}
##echo cross_subj_perf.sh GLM 21, TU, TU - trial: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
##
##sleep 1
#
#
#
