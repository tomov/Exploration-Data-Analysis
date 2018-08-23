function [X, names] = scratch

res = 100;
EXPT = exploration_expt();
hrf = spm_hrf(1 / res);

multi = exploration_create_multi(1, 1, 1);
ons = multi.onsets{1} * res;

x = zeros(round(max(ons)),1);
x(round(ons)) = 1;



[X, names] = ccnl_get_design();

figure;
plot_idx = 1;

    subplot(nplots,1,plot_idx);
    plot_idx = plot_idx + 1;
    plot(x);


function [X, names] = get_design(multi)
    % find max onset
    max_onset = 0;
    for j = 1:length(multi.names)
        max_onset = max(max_onset, max(multi.onsets{j}));
    end
    max_onset = max_onset * res;

    % iterate over events
    names = {};
    for j = 1:length(multi.names)
        onsets = multi.onsets{j} * res;

        x = zeros(ceil(max_onset),1);
        x(round(onsets)) = 1;
        x = convolve_and_subsample(x, res);
        if j == 1
            X = x;
        else
            X = [X x];
        end

        % iterate over pmods for event
        if isfield(multi, 'pmod') && j <= length(multi.pmod)
            for k = 1:length(multi.pmod(j).name)
                assert(length(onsets) == length(multi.pmod(j).param{k}), ['multi.pmod(' + num2str(j) + ').param{' + num2str(k) + '} has the wrong number of elements']);
                x = zeros(round(max(onsets)),1);
                x(round(onsets)) = multi.pmod(j).param{k};
                x = convolve_and_subsample(x, res);
                X = [X x];
            end
        end
    end

end



function x = convolve_and_subsample(x, res)
    hrf = spm_hrf(1 / res);
    x = conv(x, hrf);
    x = x(1:res:end);
end
