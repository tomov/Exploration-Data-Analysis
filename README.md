# Exploration_Exploitation

A model to understand exploration exploitation dilemma approaches in human decisions along with GLMs for fMRI correlation

Scripts are in `/ncf/gershman/Lab/scripts/matlab/Exploration/`.

Data are in `/ncf/gershman/Lab/Exploration/`.

Useful links:
- [CBS cluster FAQ](http://cbs.fas.harvard.edu/science/core-facilities/neuroimaging/information-investigators/faq) -- how to use the cluster, send jobs, ArcGet.py, slurm, sacct, etc
- [CBS central login](http://cbscentral.rc.fas.harvard.edu) -- where the fMRI data live

## To preprocess a newly scanned subject

0. Copy everything to the cluster and login
   * See `scp_to_ncf.sh`, maybe edit it accordingly and run it `./scp_to_ncf.sh`
1. Open `downloadConvertSBatch.sh` and edit `subjects` to include the **new subject only**.
2. Open `downloadConvert.sh` and edit `fileNames` if necessary. Read comments for details
3. Run `./downloadConvertSBatch.sh`
   * Takes ~30 mins to complete.
   * Make sure there are no scary errors in the .err file or the .out file.
   * Make sure all files (struct.nii, run001.nii, etc) are written normally in the .out file
   * Go to `/ncf/gershman/Lab/Exploration/subjects/` and make sure all subject files are there
4. Open `exploration_getSubjectsDirsAndRuns.m` and **append** new subject info to `subjects`, `subjdirs`, and `nRuns` accordingly (so they include info for **all subjects**). Read comments for details.
5. Open `ccnl_fmri_preproc.sh` and edit `subjects` to include the index of the **new subject only**. This index is the ordinal of the subject info in `exploration_getSubjectsDirsAndRuns.m`
6. Run `./ccnl_fmri_preproc.sh`
   * Takes ~10 hours (!) to complete.
   * Make sure there are no scary errors in the .err file or the .out file.
   * Make sure all files are in the data directory (compare with previous subjects).
7. Open MATLAB with a GUI and cd into scripts directory
   * Log into one of the cluster terminals in Northwest (e.g. by the fMRI scanner),
   * or download X11 on your Mac, `ssh -X` into the cluster, and run `matlab` from there (super slow)
8. Run `ccnl_plot_movement(exploration_expt(), XX)`, where XX is the subject index (e.g. 1)
   * Make sure subject didn't move too much during runs
9. Run `ccnl_check_registration(exploration_expt(), XX)`, where XX is the subject index
   * Click around (e.g. around ventricles), make sure average functional (top) is aligned with the structural (bottom)

## To create and run a new GLM

1. Add the GLM to `exploration_create_multi.m`
2. Test the GLM locally
    * Call `multi = exploration_create_multi(...)` in MATLAB
    * Inspect `multi` and make sure all the names, onsets, and pmods make sense
    * Batch test for all subjects and runs using `exploration_create_multi_sanity.m`
3. Copy `exploration_create_multi.m` to the cluster
    * See `scp_to_ncf.sh`
4. Log into the cluster and run an interactive job, e.g. with `./srun.sh`
5. Run MATLAB in the interactive job, e.g. with `./matlab.sh`
6. Test the GLM on the cluster
    * repeat same steps as testing locally
7. Dry run the GLM in MATLAB by calling `ccnl_fmri_glm(...)` to make sure SPM actually likes it; then interrupt it with Ctrl+C once it looks like it's working
    * Optionally use the `fake_glm=true` argument
8. Exit MATLAB and edit `ccnl_fmri_glm.sh`
    * Change `goodSubjects = (...)` to include all subjects
    * Change `for model in {...}` to use your GLM only
9. Run `./ccnl_fmri_glm.sh`

## To run a contrast

1. Add your contrast to `run_ccnl_fmri_con.m`
    * Copy-paste one of the calls to `ccnl_fmri_con` and modify accordingly
    * Please keep all previous contrasts, just comment them out
2. Copy `run_ccnl_fmri_con.m` to the cluster
3. Log into the cluster and run an interactive job, e.g. with `./srun.sh`
4. Run MATLAB in the interactive job, e.g. with `./matlab.sh`
5. Test the contrast on the cluster
    * Call `run_ccnl_fmri_con` in MATLAB, make sure it starts doing stuff and Ctrl+C it
3. Exit MATLAB and run `./run_ccnl_fmri_con.sh`

## To adapt pipeline to a new experiment

1. TODO
