% Hacky way to create a cvpartition object for given folds
% b/c MATLAB does not support this
% also have to keep the hacked cvpartition*.m files around
%
function c = cvpartition_from_folds(foldid)
    c = cvpartition(numel(foldid), 'Kfold', numel(unique(foldid)));
    c.Impl.indices = foldid;
    c.Impl.TestSize = accumarray(foldid, 1)';
    c.Impl.TrainSize = size(foldid, 1) - c.Impl.TestSize;
end
