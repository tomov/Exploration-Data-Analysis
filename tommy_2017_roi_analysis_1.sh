mkdir output

outfileprefix="output/tommy_2017_roi_analysis_1"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running tommy_2017_roi_analysis_1 GLM 21 >> jobs.txt
echo ---------------- >> jobs.txt

sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'tommy_2017_roi_analysis_1;exit'"`
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo tommy_2017_roi_analysis_1.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

