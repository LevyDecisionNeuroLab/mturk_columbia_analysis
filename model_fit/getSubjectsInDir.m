function [ subjects, datanames ] = getSubjectsInDir(data_path, prefix)
%GETSUBJECTSINDIR Extracts all subject IDs from a DATA_PATH,
%Data in .mat file are generated when running the task script. We save the data from the same subject in the folder named 'subjNo.' such as 'subj78'. 
%assuming that all folders follow the naming convention "prefix{ID}"

% get all files
subj_dirs = dir([fullfile(data_path, prefix) '*.csv']);

% subjects = cell(1, length(subj_dirs));
subjects = repmat({''}, 1, length(subj_dirs));
datanames = repmat({''}, 1, length(subj_dirs));

for k = 1 : length(subj_dirs)
    % Extract subject ids (MTurk worker ID) from folder names

    % example name: risk_3R5OYNIC3FVVK9L2EGZ19I06R48TP0_
    % ANKDLLQHHM2OH_3RANCT1ZVJ3B3CWBLICC8EC5RZ1BU7_results.csv
    ids = split(subj_dirs(k).name, '_');
    subjects{k} = ids{3};
    datanames{k} = subj_dirs(k).name;

end

% sort
subjects = sort(subjects);

end