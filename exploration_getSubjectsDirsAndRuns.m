function [ subjects, subjdirs, goodRuns, goodSubjects ] = exploration_getSubjectsDirsAndRuns()

% Get the list of subjects, subject directories, and number of runs for the
% fMRI GLM code
%

%[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% the participant id as entered in psychopy
subjects = {'uep001', 'uep002', 'uep003', 'uep004', 'uep005', ...
            'uep006', 'uep007', 'uep008', 'uep009', 'uep010', ... 
            'uep011', 'uep012', 'uep013', 'uep014', 'uep015', ...
            'uep016', 'uep017', 'uep018', 'uep019', 'uep020', ...
            'uep021', 'uep022', 'uep023', 'uep024', 'uep025', ...
            'uep026', 'uep027', 'uep028', 'uep029', 'uep030', ...
            'uep031'};

% should be identical to the list of subjects in the csv file
% and in the same order
% this is a basic assumption for getGoodSubjects() to work
% we are listing them here explicitly as a sanity check for the csv file
%
%assert(mean(strcmp(subjects, unique(data.participant)')) == 1);

% the names of the CORRESPONDING directories from CBS central
subjdirs = {'180725_UEP_001', '180727_UEP_002', '180727_UEP_003',  '180730_UEP_004', '180730_UEP_005', ...
            '180801_UEP_006', '180802_UEP_007', '180803_UEP_008',  '180803_UEP_009', '180804_UEP010',  ...
            '180804_UEP_011', '180804_UEP_012', '180804_UEP_013',  '180804_UEP_014', '180804_UEP_015', ...
            '180805_UEP_016', '180805_UEP_017', '180805_UEP_018_2','180805_UEP_019', '180805_UEP_020', ...
            '180805_UEP_021', '180806_UEP_022', '180806_UEP_023',  '180807_UEP_024', '180807_UEP_025', ...
            '180807_UEP_026', '180808_UEP_027', '180808_UEP_028',  '180808_UEP_029', '180809_UEP_030', ...
            '180809_UEP_031'};

% assumes runs are always in order: 1,2,3,4,...
%nRuns = {8,8}; % runs per subject

% which runs to include/exclude for each subject
goodRuns = {logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([0 1 0 1 1 0 1 1]), logical([1 1 1 1 0 1 1 0]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 0 1 1 1 0]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 0 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 0 0 1 1]), logical([1 1 0 0 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 0 1 1 1]), logical([1 0 1 0 0 1 1 0]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), logical([1 1 1 1 1 1 1 1]), ...
                logical([1 1 1 1 1 1 1 1]) };




% which subjects are good
goodSubjects = 1:20;
 
assert(numel(subjects) == numel(subjdirs));
assert(numel(subjects) == numel(goodRuns));
