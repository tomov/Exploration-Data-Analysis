mkdir output

outfileprefix="output/preload_betas"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running preload_betas  >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf_holy --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'preload_betas_from_masks({\'masks/mask.nii\'}, \'trial_onset\');exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo preload_betas.sh mask.nii trial_onset: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

