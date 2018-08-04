function exploration_create_multi_sanity(glmodels, subjs, runs)

% Sanity test context_create_multi.m before shipping to NCF by running
% through the GLMs and generating multi structures with differnet subjects
% / runs
% I would suggest also trying subjects that are "bad" (not in
% getGoodSubjects()), e.g. #9
% e.g. context_create_multi_sanity(154, getGoodSubjects(), 1:9)
%
% INPUT:
% glms = glmodels to test, e.g. 1:20
% subjs (optional) = subject indices to test, e.g. getGoodSubjects().
%                    Defaults to all good subjects
% runs (optional) = runs to test, e.g. 1:9. Defaults to 1:9
%

% set default parameters
%
if nargin < 3
    runs = 1:7;
end

if nargin < 2
    subjs = exploration_GoodSubjects();
end

for glmodel = glmodels
    for subj = subjs
        for run = runs
            multi = exploration_create_multi(glmodel, subj, run);
        end
    end
end