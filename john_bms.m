clear res;

load fit_AU_25nstarts_fixed.mat
res(1) = results;
load fit_ACU_25nstarts_fixed.mat
res(2) = results;
load fit_OpAL_25nstarts_fixed.mat
res(3) = results;

bms = mfit_bms(res)


clear res;

load fit_AU_25nstarts_mixed.mat
res(1) = results;
load fit_ACU_25nstarts_mixed.mat
res(2) = results;
load fit_OpAL_25nstarts_mixed.mat
res(3) = results;

bms = mfit_bms(res)


clear res;

load fit_AU_25nstarts_random.mat
res(1) = results;
load fit_ACU_25nstarts_random.mat
res(2) = results;
load fit_OpAL_25nstarts_random.mat
res(3) = results;

bms = mfit_bms(res)
