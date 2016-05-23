function [variables] = read_behavior_score(parameters)
    Score_data = importdata(parameters.score_file);
    ScoreIdx = -1;
    CtrlIdx = -1;
    for ni=2:size(Score_data.textdata, 2)
        if(length(strfind(Score_data.textdata{1,ni}, parameters.score_name)))
            ScoreIdx = ni - 1;
        else if(length(strfind(Score_data.textdata{1,ni}, parameters.control_variable_name)))
            CtrlIdx = ni - 1;
            end
        end    
    end
    variables.SubjectID = Score_data.textdata(2:end,1);    

    if(~isempty(parameters.control_variable_name))
        if(ScoreIdx < 0)
            error('Error: can not find \"%s\" in the given behavior score file: %s!\n', parameters.score_name, parameters.score_file);
        else
            variables.one_score = Score_data.data(:,ScoreIdx);
        end
        if(CtrlIdx < 0)
            error('Error: can not find \"%s\" in the given behavior score file: %s!\n', parameters.control_variable_name, parameters.score_file);
        else
            variables.control_var = Score_data.data(:, CtrlIdx);
        end         
        % remove subjects without behavior score
        Remove_idx = find(isnan(variables.control_var) + isnan(variables.one_score)~=0);
        
        if(length(Remove_idx)==1)
            warning('Subject \"%s\" was not included in the analysis because lack of behavior score or control variable:', variables.SubjectID(Remove_idx));
        else if(length(Remove_idx)>1)
                warning('The following subjects have been removed from the analysis because they do not have behavior score or control variables:\n', length(Remove_idx));
                fprintf('%s\n', variables.SubjectID{Remove_idx});
            end
        end
        
        variables.one_score(Remove_idx) = [];
        variables.SubjectID(Remove_idx) = [];
        variables.control_var(Remove_idx) = [];     
    else 
        if(ScoreIdx < 0)
            error('Error: can not find \"%s\" in the given behavior score file: %s!\n', parameters.score_name, parameters.score_file);
        else
            variables.one_score = Score_data.data(:,ScoreIdx);
        end
    
        % remove subjects without behavior score
        Remove_idx = find(isnan(variables.one_score)~=0);
        
        if(length(Remove_idx)==1)
            warning('Subject \"%s\" was not included in the analysis because lack of behavior score:', variables.SubjectID(Remove_idx));
        else if(length(Remove_idx)>1)
                warning('The following subjects have been removed from the analysis because they do not have behavior score \ns', length(Remove_idx));
                fprintf('%s\n', variables.SubjectID{Remove_idx});
            end
        end
        
        variables.one_score(Remove_idx) = [];
        variables.SubjectID(Remove_idx) = [];
    end


    variables.SubNum = length(variables.SubjectID); % number of subject used in the analysis

%     if(length(variables.one_score) ~= parameters.total_valid_record_number )
%         error('%d records was readed from the score file, is it correct?\n', length(variables.one_score));
%     end
end