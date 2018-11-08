mkdir output

outfileprefix="output/univariate_decoder"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running univariate_decoder  >> jobs.txt
echo ---------------- >> jobs.txt

echo -- RU -- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'tommy\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, RU, tommy, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'RU - trial\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, RU, RU - trial, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'TU - trial\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, RU, TU - trial, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1


echo -- TU -- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'badre\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, TU, badre, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'tommy\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, TU, tommy, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'RU - trial\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, TU, RU - trial, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'TU - trial\', 4, true, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo univariate_decoder.sh GLM 21, TU, TU - trial, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'RU - trial\', 4, true, 1);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, RU - trial, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'RU\', 4, true, 1);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, RU, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1

#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 0.01);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=0.01 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 0.1);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=0.1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 1);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 10);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=10 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'RU\', \'badre\', 4, true, 100);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, RU, badre, norm=4, do_orth=1, lambda=100 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1



#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'TU - trial\', 4, true, 0.01);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, TU, TU - trial, norm=4, do_orth=1, lambda=0.01 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'TU - trial\', 4, true, 0.1);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, TU, TU - trial, norm=4, do_orth=1, lambda=0.1 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'TU - trial\', 4, true, 10);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, TU, TU - trial, norm=4, do_orth=1, lambda=10 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1
#
#sbatch_output=`sbatch -p ncf --mem 50001 -t 10-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'univariate_decoder(21, \'TU\', \'TU - trial\', 4, true, 100);exit'"`
#sbatch_output_split=($sbatch_output)
#job_id=${sbatch_output_split[3]}
#echo univariate_decoder.sh GLM 21, TU, TU - trial, norm=4, do_orth=1, lambda=100 : ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
#
#sleep 1