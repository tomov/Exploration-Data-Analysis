#!/bin/bash
#SBATCH -p ncf # partition (queue)
#SBATCH --mem 8000 # memory
#SBATCH -t 0-8:00 # time (D-HH:MM)

module load mri_convert/2015_12_03-ncf

# download all structural and functional scans in order in which they were taken
#
ArcGet.py -a cbscentral -s ${1} -r MEMPRAGE\ RMS,Minn_HCP_1.5mm_Task

experiment="Exploration" # must match directory name

# fileNames should correspond to order in which scans were taken, e.g. normally it's fileNames=(struct run001 run002 run003 run004 run005 run006 run007 run008)
# HOWEVER, if things are out of order or some scans were bad, might have to reorder and include dummy entries for the bad scans,
# e.g. if the structural was bad and you took another structural in the end, it would look like this: fileNames=(struct_bad run001 run002 run003 run004 run005 run006 run007 run008 struct)
fileNames=(struct run001 run002 run003 run004 run005 run006 run007 run008)


cd /ncf/gershman/Lab/${experiment}/subjects/

mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/
mkdir  /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW
mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/preproc

#ArcGet.py -a cbscentral -s ${1} 

cd /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW

myruns=`ls *.MR.Investigators_Gershman.*.1.* | sort -t. -nk4`

count=0
for run in ${myruns[*]}; do
	echo ------------------
	echo Run = ${run}
	echo count = ${count}
	echo filename = ${fileNames[$count]}
	mri_convert -it siemens_dicom -ot nii -i ${run} -o /ncf/gershman/Lab/${experiment}/subjects/${1}/preproc/${fileNames[$count]}.nii
	count=${count}+1
done
