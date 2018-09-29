% behavioral BMS for AU, ACU, OpAL
%

clear res;

load fit_AU_25nstarts_fixed.mat
res(1) = results;
load fit_ACU_25nstarts_fixed.mat
res(2) = results;
load fit_OpAL_25nstarts_fixed.mat
res(3) = results;

disp('fixed effects');

%disp('AU');
%array2table(res(1).x, 'VariableNames', {res(1).param.name})
%disp('ACU');
%array2table(res(2).x, 'VariableNames', {res(2).param.name})
%disp('OpAL');
%array2table(res(3).x, 'VariableNames', {res(3).param.name})

bms = mfit_bms(res)

clear res;

load fit_AU_25nstarts_mixed.mat
res(1) = results;
load fit_ACU_25nstarts_mixed.mat
res(2) = results;
load fit_OpAL_25nstarts_mixed.mat
res(3) = results;

disp('mixed effects');

%disp('AU');
%array2table(res(1).x, 'VariableNames', {res(1).param.name})
%disp('ACU');
%array2table(res(2).x, 'VariableNames', {res(2).param.name})
%disp('OpAL');
%array2table(res(3).x, 'VariableNames', {res(3).param.name})

bms = mfit_bms(res)


clear res;

load fit_AU_25nstarts_random.mat
res(1) = results;
load fit_ACU_25nstarts_random.mat
res(2) = results;
load fit_OpAL_25nstarts_random.mat
res(3) = results;

disp('random effects');

%disp('AU');
%array2table(res(1).x, 'VariableNames', {res(1).param.name})
%disp('ACU');
%array2table(res(2).x, 'VariableNames', {res(2).param.name})
%disp('OpAL');
%array2table(res(3).x, 'VariableNames', {res(3).param.name})

bms = mfit_bms(res)


clear res;
load trash/fits_G0_N0/fit_AU_25nstarts_fixed.mat
res(1) = results;
load fit_AU_25nstarts_fixed.mat
res(2) = results;
disp('AU fixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_ACU_25nstarts_fixed.mat
res(1) = results;
load fit_ACU_25nstarts_fixed.mat
res(2) = results;
disp('ACU fixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_OpAL_25nstarts_fixed.mat
res(1) = results;
load fit_OpAL_25nstarts_fixed.mat
res(2) = results;
disp('OpAL fixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)


clear res;
load trash/fits_G0_N0/fit_AU_25nstarts_mixed.mat
res(1) = results;
load fit_AU_25nstarts_mixed.mat
res(2) = results;
disp('AU mixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_ACU_25nstarts_mixed.mat
res(1) = results;
load fit_ACU_25nstarts_mixed.mat
res(2) = results;
disp('ACU mixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_OpAL_25nstarts_mixed.mat
res(1) = results;
load fit_OpAL_25nstarts_mixed.mat
res(2) = results;
disp('OpAL mixed: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)


clear res;
load trash/fits_G0_N0/fit_AU_25nstarts_random.mat
res(1) = results;
load fit_AU_25nstarts_random.mat
res(2) = results;
disp('AU random: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_ACU_25nstarts_random.mat
res(1) = results;
load fit_ACU_25nstarts_random.mat
res(2) = results;
disp('ACU random: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)

clear res;
load trash/fits_G0_N0/fit_OpAL_25nstarts_random.mat
res(1) = results;
load fit_OpAL_25nstarts_random.mat
res(2) = results;
disp('OpAL random: old fits (G0, N0) vs. new fits');
bms = mfit_bms(res)
