#!/bin/bash
#SBATCH -p ncf # partition (queue)
#SBATCH --mem 8000 # memory
#SBATCH -t 0-8:00 # time (D-HH:MM)



ArcGet.py -a cbscentral -s ${1} -r MEMPRAGE\ RMS,Minn_HCP_1.5mm_Task
experiment="Exploration"

#Note that struct needs to come first
fileNames=(struct run001 run002 run003 run004 run005 run006 run007 run008) 


cd /ncf/gershman/Lab/${experiment}/subjects/

if [ -d "/ncf/gershman/Lab/${experiment}/subjects/${1}/" ];
then
	echo "directory already exists"
else
	mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/
fi

mkdir  /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW
cd /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW
mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/preproc

#ArcGet.py -a cbscentral -s ${1} 

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
