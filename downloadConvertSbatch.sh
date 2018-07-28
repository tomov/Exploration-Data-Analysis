#!/bin/bash

module load mri_convert/2015_12_03-ncf

experiment="Exploration"

subjects=('180725_UEP_001')



#subjects=$(/ncf/gershman/Lab/${experiment}/subjects.txt)

for subj in ${subjects[*]}; do

echo $subj
cd /ncf/gershman/Lab/${experiment}/subjects/
mkdir /ncf/gershman/Lab/${experiment}/subjects/${subj}/
sbatch -o /ncf/gershman/Lab/${experiment}/subjects/${subj}_%j.out -e /ncf/gershman/Lab/${experiment}/subjects/${subj}_%j.err /ncf/gershman/Lab/scripts/matlab/Exploration/downloadConvert.sh ${subj}
done

