mkdir output

outfileprefix="output/univariate_decoder_refactored"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running univariate_decoder_refactored  >> jobs.txt
echo ---------------- >> jobs.txt

#function univariate_decoder_refactored(roi_glmodel, roi_contrast, glmodel, regressor, do_orth, lambda, standardize, mixed_effects, clusterFWEcorrect, extent, Num, intercept, flip_sign, do_CV, get_null)

# GLM 36 V/TU, extent 100, mixed, flip, intercept on
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(45, \'RU\', 45, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(45, \'TU\', 45, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(46, \'RU\', 46, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(46, \'TU\', 46, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 46, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 46, \'VTU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# conjunction (avg) V/TU, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1, intercept on/off
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(-1, \'conj_3\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 45, \'VTU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 64, \'VTU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_36\', 36, \'VTU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     )

# DV
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(29, \'DV\', 29, \'DV\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# DV, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1, intercept on
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(29, \'DV\', 29, \'DV\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(29, \'DV\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(29, \'DV\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# conjunction (avg)j, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1, intercept on
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(-1, \'conj_3\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 64, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_3\', 64, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_36\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'conj_36\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# GLM 66 badre, mixed, flip_sign (as it should be), lambda = 1, intercept on/off
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(-1, \'badre\', 66, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 66, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 66, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 66, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     )

# GLM 45, badre, mixed, flip_sign (as it should be), lambda = 1, intercept on
declare -a fn_calls=(
                     "univariate_decoder_refactored(-1, \'badre\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
                     "univariate_decoder_refactored(-1, \'badre\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
                     )

# GLM 45, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1, intercept off
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(45, \'RU\', 45, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(45, \'RU\', 45, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(45, \'TU\', 45, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(45, \'TU\', 45, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     )

# separate GLMs, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1, intercept on
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(39, \'RU\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(40, \'RU\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(39, \'TU\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(40, \'TU\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# GLM 36, extent >= 100, mixed eff, flip_sign (as it should be), lambda = 1/40 (as in fitrlinear)
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', false, 1/40, 2, true, false, 100, 1, false, true, false, false)"
#                     )



#declare -a fn_calls=(
#                     "univariate_decoder_refactored(-1, \'badre\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# single-regressor GLMs
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(39, \'RU\', 39, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(40, \'TU\', 40, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )

# cluster FWE, Num = 3
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, true, true, [], 3, true, true, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1, 2, true, true, [], 3, true, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1, 2, true, true, [], 3, true, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, true, true, [], 3, true, true, false, false)"
#                     )

# extent >= 100, fixed eff, single-regressor
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(39, \'RU\', 39, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(40, \'TU\', 40, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 39, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 40, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     )

# extent >= 100, mixed eff, flipped sign, NO INTERCEPT
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, true, false, false)"
#                     )

# extent >= 100, mixed eff, NO flipped sign, NO INTERCEPT
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', false, 1, 2, true, false, 100, 1, false, false, false, false)"
#                     )


# extent >= 100, fixed eff, flipped sign
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, true, false, false)"
#                     )

# extent >= 100, fixed eff (as in paper)
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', false, 1, 2, false, false, 100, 1, false, false, false, false)"
#                     )





# controls, extent >= 100, mixed eff
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(64, \'RU\', 64, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(64, \'RU\', 64, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(64, \'TU\', 64, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(64, \'TU\', 64, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 64, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 64, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 64, \'RU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 64, \'TU\', false, 1, 2, true, false, 100, 1, true, true, false, false)"
#                     )
#

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, true)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 0.0001, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 0.001, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 0.01, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 0.1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 10, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 100, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 1000, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 10000, 2, false, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(47, \'DV\', 47, \'DV\', false, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )



#declare -a fn_calls=(
#                     "univariate_decoder_refactored(45, \'RU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(45, \'RU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(45, \'TU\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(45, \'TU\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 45, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1, false)"
#                     )

# repro from paper
#declare -a fn_calls=(
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'badre\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'RU\', true, 1, 2, false, false, 100, 1, false)"
#                     "univariate_decoder_refactored(-1, \'dlpfc\', 36, \'TU\', true, 1, 2, false, false, 100, 1, false)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 10, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 100, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 1000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', false, 10000, 2, true, false, 100, 1, true)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', false, 10000, 2, true, false, 100, 1, true)"
#                     )


#declare -a fn_calls=(
#                     "univariate_decoder_refactored(41, \'V\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(41, \'V\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(41, \'V\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'RU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'RU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'TU\', 41, \'V\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'RU\', true, 1, 2, true, false, 100, 1)"
#                     "univariate_decoder_refactored(36, \'TU\', 36, \'TU\', true, 1, 2, true, false, 100, 1)"
#                     )

#declare -a fn_calls=(
#                     "univariate_decoder_refactored(11, \'RU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder_refactored(11, \'RU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder_refactored(11, \'TU\', 35, \'RU\', true, 1, 2, false, false, 100)"
#                     "univariate_decoder_refactored(11, \'TU\', 35, \'TU\', true, 1, 2, false, false, 100)"
#                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo univariate_decoder_refactored.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

