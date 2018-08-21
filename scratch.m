res = 100;

multi = exploration_create_multi(1, 1, 1);
ons = multi.onsets{1} * res;

x = zeros(round(max(ons)),1);
x(round(ons)) = 1;

figure;

subplot(3,1,1);
plot(x(1:EXPT.TR*res:end));

subplot(3,1,2);
plot(hrf(1:EXPT.TR*res:end));

subplot(3,1,3);
y = conv(x,hrf);
plot(y(1:EXPT.TR*res:end));

max_onset = 0;
for j = 1:length(multi.names)
    max_onset = max(max_onset, max(multi.onsets{j}));
end
max_onset = max_onset * res;

names = {};
for j = 1:length(multi.names)
    onsets = multi.onsets{j} * res;

    x = zeros(ceil(max_onset),1);
    x(round(onsets)) = 1;
    x = convolve_and_subsample(x);
    if j == 1
        X = x(1:EXPT.TR*res:end);
    else
        X = [X x(1:EXPT.TR*res:end)];
    end

    if isfield(multi, 'pmod') && j <= length(multi.pmod)
        for k = 1:length(multi.pmod(j).name)
            assert(length(onsets) == length(multi.pmod(j).param{k}), ['multi.pmod(' + num2str(j) + ').param{' + num2str(k) + '} has the wrong number of elements']);
            x = zeros(round(max(onsets)),1);
            x(round(onsets)) = multi.pmod(j).param{k};
            X = [X x(1:EXPT.TR*res:end)];
        end
    end
end


function x = convolve_and_subsample(x, TR, res)
    hrf = spm_hrf(1 / res);
    x = conv(x, hrf);
    x = x(1:TR*res:end);
end
