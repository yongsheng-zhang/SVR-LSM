function [output_score] = regressout_variables(input_score, control_variable)
    
    tmp_size = size(control_variable);
    if(tmp_size(1)>tmp_size(2))
        control_variable = control_variable.';
    end
    
    tmp_size = size(input_score);
    if(tmp_size(1)>tmp_size(2))
        input_score = input_score.';
    end
    
    covy12=cov(input_score,control_variable);
    C=inv([var(input_score) covy12(2,1);  covy12(2,1) var(control_variable)])*[ covy12(2,1); covy12(2,1)];
    M=C'*[input_score;control_variable];
    output_score=(input_score-M).';
end
