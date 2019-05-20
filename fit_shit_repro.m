% reproduce circularity in multilinear_fit
% 
% ... o m g..........
% I figured it out....................
% it's not circular, but of course predictions are anticorrelated across folds
% b/c there averages of the folds are significantly different from each other....
% e.g. fold 1 is very low, training on folds 2..K will learn to predict a higher average (b/c the other folds are higher) => it will make a higher prediction on fold 1
% convsersely, if fold 2 is very high, training on the other folds will predict a higher average => it will make a lower prediction on fold 2
% and so on...
%
% => cannot use correlation across folds to assess fit
% maybe use correlation within folds?

load fit_shit.mat




cv = cvpartition_from_folds(foldid);

kfold = length(unique(foldid));
for k = 1:kfold
    which = training(cv, k);
    pred(~which,:) = ridgepred(X(which,:), y(which,:), X(~which,:), Lambda);
    yavg(~which,:) = mean(y(~which));
end
mse = immse(pred, y);

[r,p] = corr(y, pred)


plot(pred);
hold on;
plot(-yavg);
hold off;
legend({'pred', '-yavg'});
