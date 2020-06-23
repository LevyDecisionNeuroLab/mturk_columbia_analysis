% NOTE: Requires MATLAB optim library
% notice to change the constrained/unconstrained function, and change the fitpat.mat file name to constrained/unconstrained

clearvars
close all

%% Input set up
fitparwave = 'data_test_3sub';
search = 'grid'; % which method for searching optimal parameters
model = 'ambigNrisk'; % which utility function
isconstrained = 0; % if use constrained fitting. 0-unconstrained, 1-constrained, 2-both
includeAmbig = 1; % whether to include ambiguous trials or not

%% Set up loading + subject selection
% TODO: Maybe grab & save condition somewhere?

root = 'Z:\Lab_Projects\mturk_Columbia\behavioral\'; % Need to change if doing analysis in different folders
data_path = fullfile(root, fitparwave); % root of folders is sufficient
fitpar_out_path = fullfile(root,'model_fit_results', fitparwave);
%graph_out_path  = fullfile(root, 'ChoiceGraphs/');

summary_file = fullfile(fitpar_out_path, ['param_nonparam_' fitparwave '.txt']); % parametric and nonparametric risk and ambig attitudes
choiceData_file = fullfile(fitpar_out_path, ['choice_data_' fitparwave '.xls']); % choice matrix

if exist(fitpar_out_path)==0
    mkdir(fullfile(root,'model_fit_results'),fitparwave)
end

addpath(genpath(data_path)); % generate path for all the subject data folder

% get MTurker IDs and file names
[subjects, datanames] = getSubjectsInDir(data_path, 'risk');

% exclude subjects
% exclude = []; 
% subjects = subjects(~ismember(subjects, exclude));

% for refitting the subjects needing constraints

% poolobj = parpool('local', 8);

%%
tic

for subj_idx = 1:length(subjects)
    %% read data and clean   
    subjectNum = subjects{subj_idx};
    
    dataname = datanames{subj_idx};
    
    Data = readtable(dataname);

    choice_data = Data(1:height(Data)-1, 1:8);
    bonus_data = Data(end, 9:13); 

    %% Refine variables
   
    % trial paramters
    values_all = choice_data.trial_reward;
    ambigs_all = choice_data.trial_uncertainty;
    probs_all  = choice_data.trial_uncertainty;
    % correct uncertainty
    ambigs_all(strcmp(choice_data.trial_type, 'risk')) = 0;
    probs_all(strcmp(choice_data.trial_type, 'ambiguity')) = 0.5;

    % choices
    % !TODO Check if there are missing responses
    choices_all = zeros(size(ambigs_all));    
    choices_all(strcmp(choice_data.user_choice, 'risk'))=1;
           
    %% Clean data 
    
    % get rid of trials with reward less than $5, and missing responses
    trial_mask = values_all>=5 & ~isnan(choices_all);
    
    choice = choices_all(trial_mask);
    values = values_all(trial_mask);
    ambigs = ambigs_all(trial_mask);
    probs = probs_all(trial_mask);
    
    %% Prepare variables for model fitting

    fixed_valueP = 5; % Value of fixed reward
    fixed_prob = 1;   % prb of fixed reward 
    ambig = unique(ambigs(ambigs > 0)); % All non-zero ambiguity levels 
    prob = unique(probs); % All probability levels
    base = 0; % another parm in the model. Not used.

    if strcmp(search, 'grid')
    % grid search
    % range of each parameter
        if strcmp(model,'ambigNrisk')
            slopeRange = -4:0.2:1;
            bRange = -2:0.2:2;
            aRange = 0:0.2:4;
        else
            slopeRange = -4:0.2:1;
            bRange = -2:0.2:2;
            aRange = -2:0.2:2;
        end
        % three dimenstions
        [b1, b2, b3] = ndgrid(slopeRange, bRange, aRange);
        % all posibile combinatinos of three parameters
        b0 = [b1(:) b2(:) b3(:)];
    elseif strcmp(search,'single')
        % single search
        b0 = [-1 0.5 0.5]; % starting point of the search process, [gamma, beta, alpha]
    elseif strcmp(search, 'random')
        % independently randomized multiple search starting points
        bstart = [-1 0 1]; % starting point of the search process, [gamma, beta, alpha]
        itr = 100; % 100 iteration of starting point
        b0 = zeros(itr,length(bstart));
        for i = 1:itr
            % gamma: negative, around -1, so (-2,0)
            % beta: [-1,1] possible to be larger than 1?
            % alpha: (0,4)
            b0(i,:) = bstart + [-1+2*rand(1) -1+2*rand(1) -1+2*rand(1)]; % randomize search starting point, slope, beta, alpha
        end
    end


    refVal = fixed_valueP * ones(length(choice), 1);
    refProb = fixed_prob  * ones(length(choice), 1);        

    %% Fit model

    % Two versions of function, calculate both the unconstrained and constrained fittings:
    % fit_ambgiNrisk_model: unconstrained
    if isconstrained == 0 || isconstrained == 2
        [info_uncstr, p_uncstr] = fit_ambigNrisk_model(choice', ...
            refVal', ...
            values', ...
            refProb', ...
            probs', ...
            ambigs', ...
            model, ...
            b0, ...
            base);

        slope_uncstr = info_uncstr.b(1);
        a_uncstr = info_uncstr.b(3);
        b_uncstr = info_uncstr.b(2);
        r2_uncstr = info_uncstr.r2;

        disp(['Subject ' subjectNum ' unconstrained fitting completed'])

    end

    if isconstrained == 1 || isconstrained == 2
        % fit_ambigNrisk_model_Constrained: constrained on alpha and beta    
        [info_cstr, p_cstr] = fit_ambigNrisk_model_Constrained(choice', ...
            refVal', ...
            values', ...
            refProb', ...
            probs', ...
            ambigs', ...
            model, ...
            b0, ...
            base);

        slope_cstr = info_cstr.b(1);
        a_cstr = info_cstr.b(3);
        b_cstr = info_cstr.b(2);
        r2_cstr = info_cstr.r2;

        disp(['Subject ' subjectNum ' constrained fitting completed'])
    end

    %% Create choice matrices

    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.

    % Inputs: 
    %  Data
    %   .values, .ambigs, .probs, .choices (filtered by include_indices and transformed)
    %  ambig, prob (which are subsets of ambigs and probs, ran through `unique`)
    %
    % Outputs:
    %  ambigChoicesP
    %  riskyChoicesP
    %
    % Side-effects:
    %  one graph generated per-subject-domain
    %  .ambigChoicesP and .riskyChoicesP saved into `fitpar` file

    % Ambiguity levels by payoff values
    valueP = unique(values_all(ambigs_all > 0 & values_all ~= 4)); % each lottery payoff value under ambiguity
    ambigChoicesP = zeros(length(ambig), length(valueP)); % each row an ambiguity level
    for i = 1:length(ambig)
        for j = 1:length(valueP)
            selection = find(ambigs == ambig(i) & values == valueP(j));
            if ~isempty(selection)
                ambigChoicesP(i, j) = choice(selection);
            else
                ambigChoicesP(i, j) = NaN;
            end
        end
    end

    % Create riskyChoicesP
    % Risk levels by payoff values
    valueP = unique(values_all(ambigs_all == 0 & values_all ~= 4));
    riskyChoicesP = zeros(length(prob), length(valueP));
    for i = 1:length(prob)
        for j = 1:length(valueP)
            selection = find(probs == prob(i) & values == valueP(j) & ambigs == 0);
            if ~isempty(selection)
                riskyChoicesP(i, j) = choice(selection);
            else
                riskyChoicesP(i, j) = NaN;
            end
        end
    end

    riskyChoices_byLevel = zeros(1, length(prob));
    ambigChoices_byLevel = zeros(1, length(ambig));
    % Creat risky/ambig choiecs by level (nonparametric), excluding the value 5
    for i=1:length(prob)
        riskyChoices_byLevel(1,i) = nanmean(riskyChoicesP(i,2:length(riskyChoicesP)));
    end
    for i=1:length(ambig)
        ambigChoices_byLevel(1,i) = nanmean(ambigChoicesP(i,2:length(ambigChoicesP)));
    end        

    %% Write results into files
 
    % both constrained and unconstrained
    if isconstrained == 2
        % results file
%         fid = fopen([path summary_file],'w')
%         fprintf(fid,'\tPar_unconstrained\t\t\t\t\t\t\t\t\t\t\t\t\t\tPar_constrained\t\t\t\t\t\t\t\t\t\t\t\t\t\tNonPar\n')
%         fprintf(fid,'subject\tgains\t\t\t\t\t\t\tlosses\t\t\t\t\t\t\tgains\t\t\t\t\t\t\tlosses\t\t\t\t\t\t\tgains\t\t\t\t\t\tlosses\n')
%         fprintf(fid,'\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\tG_risk25\tG_risk50\tG_risk75\tG_amb24\tG_amb50\tG_amb74\tL_risk25\tL_risk50\tL_risk75\tL_amb24\tL_amb50\tL_amb74\n')
    end

    % unconstrained
    if isconstrained == 0
        % results file
        if subj_idx == 1 % write header
            fid = fopen(summary_file,'w')
            fprintf(fid,'mturk_id\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\texitflag\tmodel\toptimizer\trisk25\trisk50\trisk75\tamb24\tamb50\tamb74\n');
        end
        
        %write into param text file
        fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n',...
           subjectNum, info_uncstr.b(3), info_uncstr.b(2), info_uncstr.b(1),...
            info_uncstr.r2, info_uncstr.LL, info_uncstr.AIC, info_uncstr.BIC, ...
            info_uncstr.exitflag, info_uncstr.model, info_uncstr.optimizer,...
            riskyChoices_byLevel,ambigChoices_byLevel);

    end

    % constrained
    if isconstrained == 1
        % results file
%         fid = fopen([path summary_file],'w')
%         fprintf(fid,'\tPar_constrained\n')
%         fprintf(fid,'subject\tgains\t\t\t\t\t\t\tlosses\t\t\t\t\t\t\tgains\t\t\t\t\t\tlosses\n')
%         fprintf(fid,'\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\talpha\tbeta\tgamma\tr2\tLL\tAIC\tBIC\tG_risk25\tG_risk50\tG_risk75\tG_amb24\tG_amb50\tG_amb74\tL_risk25\tL_risk50\tL_risk75\tL_amb24\tL_amb50\tL_amb74\n')
    end
    
    choices_allP = [riskyChoicesP; ambigChoicesP];
    all_data_subject = [valueP'; choices_allP];

    dlmwrite(choiceData_file, subjectNum , '-append', 'roffset', 1, 'delimiter', ' ');  
    dlmwrite(choiceData_file, all_data_subject, 'coffset', 1, '-append', 'delimiter', '\t');
    
%     save_mat(Data, subjectNum, domain, fitpar_out_path);
end

fclose(fid)

toc

% delete(poolobj)