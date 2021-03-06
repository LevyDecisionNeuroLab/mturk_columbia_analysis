% this script calculate and print out subjective value of each lottery for
% each subjects

clearvars
close all

fitparwave = 'Behavior data fitpar_08120119';
outputwave = '_03270118';
isconstrained = 2;
% exclude should match those in the fit_parameters.m script
exclude = [77 1218]; 
% TEMPORARY: subjects incomplete data (that the script is not ready for)

%% folder and subjects
root = 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Analysis Ruonan';
data_path = fullfile(root, 'Behavior data of PTB log/'); % Original log from PTB
subjects = getSubjectsInDir(data_path, 'subj'); %function
subjects = subjects(~ismember(subjects, exclude));

path_out = fullfile(root, 'Fitpar files', fitparwave, filesep);
% cd 'D:\Ruonan\Projects in the lab\VA_RA_PTB\Analysis Ruonan\Behavior data fitpar_091017';

% defining monetary values
valueP = [4 5 6 7 8 10 12 14 16 19 23 27 31 37 44 52 61 73 86 101 120];
value = repmat(valueP,6,1);
valueN = [-4 -5 -6 -7 -8 -10 -12 -14 -16 -19 -23 -27 -31 -37 -44 -52 -61 -73 -86 -101 -120];
% six risk and ambig levels
probs = [0.25; 0.5; 0.75; 0.5; 0.5; 0.5];
ambigs = [0; 0; 0; 0.24; 0.5; 0.74];

for s = 1:length(subjects)
    subject = subjects(s);
%     subject = 53;

    % matrix for subjective value, row is riskambig levels, column is reward magnitude
    
    % constrained fitting
    if isconstrained ==1 || isconstrained == 2
        load([path_out 'RA_GAINS_' num2str(subject) '_fitpar.mat']);
        sv_constr_gains = ambig_utility(0, ...
              value, ...
              probs, ...
              ambigs, ...
              Data.alpha_cstr, ...
              Data.beta_cstr, ...
              'ambigNrisk');

        svRef_constr_gains = ambig_utility(0, 5, 1, 0, Data.alpha_cstr, Data.beta_cstr, 'ambigNrisk');

        load([path_out 'RA_LOSS_' num2str(subject) '_fitpar.mat']);
        sv_constr_loss = ambig_utility(0, ...
              value, ...
              probs, ...
              ambigs, ...
              Data.alpha_cstr, ...
              Data.beta_cstr, ...
              'ambigNrisk');

        sv_constr_loss = -1 .* sv_constr_loss;

        svRef_constr_loss = ambig_utility(0, 5, 1, 0, Data.alpha_cstr, Data.beta_cstr, 'ambigNrisk');
        svRef_constr_loss = -1 .* svRef_constr_loss;
    end
    
    % unconstrained fitting
    if isconstrained ==0 || isconstrained == 2
        load([path_out 'RA_GAINS_' num2str(subject) '_fitpar.mat']);
        sv_unconstr_gains = ambig_utility(0, ...
              value, ...
              probs, ...
              ambigs, ...
              Data.alpha_uncstr, ...
              Data.beta_uncstr, ...
              'ambigNrisk');

        svRef_unconstr_gains = ambig_utility(0, 5, 1, 0, Data.alpha_uncstr, Data.beta_uncstr, 'ambigNrisk');
        
        load([path_out 'RA_LOSS_' num2str(subject) '_fitpar.mat']);
        sv_unconstr_loss = ambig_utility(0, ...
              value, ...
              probs, ...
              ambigs, ...
              Data.alpha_uncstr, ...
              Data.beta_uncstr, ...
              'ambigNrisk');

        sv_unconstr_loss = -1 .* sv_unconstr_loss;

        svRef_unconstr_loss = ambig_utility(0, 5, 1, 0, Data.alpha_uncstr, Data.beta_uncstr, 'ambigNrisk');
        svRef_unconstr_loss = -1 .* svRef_unconstr_loss;
    end
    
    %% for Excel file - subjective values
    if isconstrained == 1 || isconstrained == 2
        xlFile = [path_out 'SV_constrained_by_lottery.xls'];
        dlmwrite(xlFile, subject, '-append', 'roffset', 1, 'delimiter', ' '); 
        dlmwrite(xlFile, svRef_constr_gains, '-append', 'coffset', 1, 'delimiter', '\t');
        dlmwrite(xlFile, sv_constr_gains, 'coffset', 1, '-append', 'delimiter', '\t');
        dlmwrite(xlFile, svRef_constr_loss, '-append', 'coffset', 1, 'delimiter', '\t');
        dlmwrite(xlFile, sv_constr_loss, 'coffset', 1, '-append', 'delimiter', '\t');
    end
    
    if isconstrained == 0 || isconstrained == 2
        xlFile = [path_out 'SV_unconstrained_by_lottery.xls'];
        dlmwrite(xlFile, subject, '-append', 'roffset', 1, 'delimiter', ' '); 
        dlmwrite(xlFile, svRef_unconstr_gains, '-append', 'coffset', 1, 'delimiter', '\t');
        dlmwrite(xlFile, sv_unconstr_gains, 'coffset', 1, '-append', 'delimiter', '\t');
        dlmwrite(xlFile, svRef_unconstr_loss, '-append', 'coffset', 1, 'delimiter', '\t');
        dlmwrite(xlFile, sv_unconstr_loss, 'coffset', 1, '-append', 'delimiter', '\t');
    end

end

