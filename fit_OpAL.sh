mkdir output

outfileprefix="output/fit_OpAL"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running fit_OpAL  >> jobs.txt
echo ---------------- >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_OpAL(load_data, 25, \'params_OpAL_25nstarts_random.mat\', 0, 0);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo fit_OpAL.sh random: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_OpAL(load_data, 25, \'params_OpAL_25nstarts_random_hierarchical.mat\', 1, 0);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo fit_OpAL.sh random hierarchical: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

# send the job to NCF
#
sbatch_output=`sbatch -p ncf --mem 50000 -t 20-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_OpAL(load_data, 25, \'params_OpAL_25nstarts_fixed.mat\', 0, 1);exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo fit_OpAL.sh fixed: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt


sleep 1