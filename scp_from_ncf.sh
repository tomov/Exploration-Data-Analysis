# Copy group-level contrasts for given model(s) from cluster to local directory.
#
# Edit the for loop (loops over models) and change the password where it says password.
# Creates a directory glmOutput in the current directory and copies stuff there.
# Make sure to install sshpass first (see below).
# WARNING: Make sure NOT to upload this to github after putting in an actual password
#

# STEP 1: install sshpass
# do it the hard way -- download from http://sourceforge.net/projects/sshpass/files/sshpass/1.05/sshpass-1.05.tar.gz
# then ./configure
# then sudo make install
# which sshpass <-- to find out where it is, and use that path in the --rsh thing below


rsync -avh --rsh="/usr/local/bin/sshpass -p password ssh -o StrictHostKeyChecking=no -l mtomov13" mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/Exploration/glmOutput/mean.nii glmOutput/mean.nii
for i in {29..29}
do
    rsync -avh --rsh="/usr/local/bin/sshpass -p password -o StrictHostKeyChecking=no -l mtomov13" mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/Exploration/glmOutput/model$i/con* glmOutput/model$i/

# Optionally copy the subject betas too (WARNING: takes up lots of memory; I recommend using an external HDD
#
#    for j in {1..31}
#    do
#        rsync -avh --rsh="/usr/local/bin/sshpass -p password ssh -o StrictHostKeyChecking=no -l mtomov13" mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/ConLearn/glmOutput/model$i/subj$j/SPM.mat glmOutput/model$i/subj$j/
#   #     rsync -avh --rsh="/usr/local/bin/sshpass -p password ssh -o StrictHostKeyChecking=no -l mtomov13" mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/ConLearn/glmOutput/model$i/subj$j/beta* glmOutput/model$i/subj$j/
#    done
done
