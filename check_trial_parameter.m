%data_path = 'Y:/mturk_Columbia/Raw_Data/Sandbox_Runs';
data_path = 'C:\Users\yl2268\Downloads\Re__Risk_Ambiguity_DONE_'
cd(data_path)

print_on_screen = true;

filename = 'check_trial_parameter_new.txt';

% get all subjects data
files = dir('risk_*_results.csv');


for sub_idx = 1:length(files)
%    sub_idx = 2;
    % read data into a table
     data_sub = readtable(files(sub_idx).name);
%    data_sub = readtable('risk_3D17ECOUOFJNU3WQES53XMTWQ1Z313_A38HVPKCR10RG1_3NG53N1RLX7FZCJS1DGRDJFLCKR8P3_results.csv');
    % clean table
    choice_data = data_sub(1:height(data_sub)-1, 1:8);
    bonus_data = data_sub(end, 9:13); % from YTL: (end, 9:13) was (121, 9:13)
    
    % convert strings to numeric
    for i = 1:height(choice_data)
        if strcmp(choice_data.trial_type{i}, 'ambiguity')
            choice_data.ambig(i) = choice_data.trial_uncertainty(i);
            choice_data.prob(i) = 0.5;
        elseif strcmp(choice_data.trial_type{i}, 'risk')
            choice_data.ambig(i) = 0;
            choice_data.prob(i) = choice_data.trial_uncertainty(i);
        end
        
        if strcmp(choice_data.user_choice{i}, 'risk') 
            choice_data.choice(i) = 1;
        elseif strcmp(choice_data.user_choice{i}, 'certainty') 
            choice_data.choice(i) = 0;
        end
    end
    
    % number of trials
    trial_num = height(choice_data);
    trial_num_ambig = sum(choice_data.ambig >0);
    trial_num_risk = sum(choice_data.ambig == 0);
    
    % ambiguity levels
    ambig_level = unique(choice_data.ambig(choice_data.ambig > 0));
    
    % risk levels
    prob_level = unique(choice_data.prob(choice_data.ambig == 0));
    
    % value levels
    value_level = unique(choice_data.trial_reward);
    
    % choice by trial table
    choice_sum_table = zeros(length(ambig_level) + length(prob_level),...
        length(value_level));
    
    repeat_table = zeros(length(ambig_level) + length(prob_level),...
        length(value_level));
    
    winning_color_table = NaN(length(ambig_level) + length(prob_level),...
        length(value_level));
    
    for trial_idx = 1:trial_num
        
        if choice_data.ambig(trial_idx) == 0 
            uncertain_idx = [prob_level; 0; 0; 0] == choice_data.prob(trial_idx);
        elseif choice_data.ambig(trial_idx) > 0
            uncertain_idx = [0; 0; 0; ambig_level] == choice_data.ambig(trial_idx);
        end
        
        value_idx = value_level == choice_data.trial_reward(trial_idx);
        
        if strcmp(choice_data.trial_winning_color{trial_idx}(1),'r')
            winning_color_table(uncertain_idx, value_idx) = 1;
        elseif strcmp(choice_data.trial_winning_color{trial_idx}(1),'b')
            winning_color_table(uncertain_idx, value_idx) = 2;
        end
        
        choice_sum_table(uncertain_idx, value_idx) = choice_sum_table(uncertain_idx, value_idx) + choice_data.choice(trial_idx);
        repeat_table(uncertain_idx, value_idx) = repeat_table(uncertain_idx, value_idx)+1;
    end
    
    choice_prob_table = choice_sum_table ./ repeat_table;
    
    % correct sum of choosing lotteries for 
    choice_sum_table(repeat_table == 0) = NaN;
    
    if print_on_screen
        % print on command window
        disp(['Subject ', num2str(sub_idx)])
        disp(files(sub_idx).name)
        disp(['No. of trials: ', num2str(trial_num)])
        disp(['No. of risky trials: ', num2str(trial_num_risk)])
        disp(['No. of ambiguous trials: ', num2str(trial_num_ambig)])
        disp('Risk levels:')
        disp(num2str(prob_level'))
        disp('Ambiguity levels:')
        disp(num2str(ambig_level'))
        disp('Reward levels:')
        disp(num2str(value_level'))
        disp('No. of reward levels')
        disp(num2str(length(value_level)))
        disp('For tables printed below, 6 rows (3 risk + 3 ambiguity levels) by 20 columns(reward levels)')
        disp('Winning color of the lottery (1red, 2blue)')
        disp(num2str(winning_color_table))
        disp('Repeats of trials:')
        disp(num2str(repeat_table))
        disp('Choice probability of choosing the lottery:')
        disp(num2str(choice_prob_table))    
        disp('Total number of times of choosing the lottery:')
        disp(num2str(choice_sum_table))        
        disp(' ')
    else        
        % print to file
        dlmwrite(filename,['Subject ', num2str(sub_idx)],'-append',...
            'delimiter','','roffset',1)
        dlmwrite(filename,['Subject ', files(sub_idx).name],'-append',...
            'delimiter','') 
        dlmwrite(filename,['No. of trials: ', num2str(trial_num)],'-append',...
            'delimiter','')
        dlmwrite(filename,['No. of risky trials: ', num2str(trial_num_risk)],'-append',...
            'delimiter','')
        dlmwrite(filename,['No. of ambiguous trials: ', num2str(trial_num_ambig)],'-append',...
            'delimiter','')   
        dlmwrite(filename,'Risk levels:','-append',...
            'delimiter','') 
        dlmwrite(filename, num2str(prob_level'),'-append',...
            'delimiter','')
        dlmwrite(filename, 'Ambiguity levels:','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(ambig_level'),'-append',...
            'delimiter','')
        dlmwrite(filename, 'Reward levels:','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(value_level'),'-append',...
            'delimiter','')
        dlmwrite(filename, 'No. of reward levels','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(length(value_level)),'-append',...
            'delimiter','')
        dlmwrite(filename, 'For tables printed below, 6 rows (3 risk + 3 ambiguity levels) by 20 columns(reward levels)','-append',...
            'delimiter','')
        dlmwrite(filename, 'Winning color of the lottery (1red, 2blue)','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(winning_color_table),'-append',...
            'delimiter','')
        dlmwrite(filename, 'Repeats of trials:','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(repeat_table),'-append',...
            'delimiter','')
        dlmwrite(filename, 'Choice probability of choosing the lottery:','-append',...
            'delimiter','')
        dlmwrite(filename, num2str(choice_prob_table),'-append',...
            'delimiter','')
        dlmwrite(filename, 'Total number of times of choosing the lottery:','-append',...
            'delimiter','')  
        dlmwrite(filename, num2str(choice_sum_table),'-append',...
            'delimiter','')
    end
    
end