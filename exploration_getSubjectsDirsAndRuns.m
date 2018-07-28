function [ subjects, subjdirs, nRuns ] = context_getSubjectsDirsAndRuns()

% Get the list of subjects, subject directories, and number of runs for the
% fMRI GLM code
%

%[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% the participant id as entered in psychopy
subjects = {'uep001', 'uep002'}; 
        
% should be identical to the list of subjects in the csv file
% and in the same order
% this is a basic assumption for getGoodSubjects() to work
% we are listing them here explicitly as a sanity check for the csv file
%
%assert(mean(strcmp(subjects, unique(data.participant)')) == 1);

% the names of the CORRESPONDING directories from CBS central
subjdirs = {'180725_UEP_001','180727_UEP_002'};

% assumes runs are always in order: 1,2,3,4,...
nRuns = {8,8}; % runs per subject



