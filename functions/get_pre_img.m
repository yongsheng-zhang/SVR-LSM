function [output, variables] = get_pre_img(parameters, variables, cmd)

    fprintf('\nBegin to run SVR-LSM...\n');
    if(~exist(['./', parameters.timestamp, '/', parameters.score_name, '/SVR_LSM']))
        eval(['mkdir ./', parameters.timestamp, '/', parameters.score_name, '/SVR_LSM']);
    end
    
    nx = variables.vo.dim(1);
    ny = variables.vo.dim(2);
    nz = variables.vo.dim(3);
    
    if(size(variables.one_score,1)<size(variables.one_score,2))
        variables.one_score = variables.one_score.';
    end
    m = svmtrain(variables.one_score,sparse(variables.lesion_dat),cmd);
    
    %% compute the beta-map
    w = m.sv_coef.'*m.SVs;
    variables.img_scale = 10/max(abs(w));
    output.beta_map = zeros(nx, ny, nz);
    output.beta_map(variables.l_idx) = w.'*variables.img_scale;
    
    variables.vo.fname = ['./', parameters.timestamp, '/', parameters.score_name, '/SVR_LSM/beta_map', parameters.surfix,'.nii'];
    spm_write_vol(variables.vo, output.beta_map);
    variables.vo.fname = ['./', parameters.timestamp, '/', parameters.score_name, '/SVR_LSM/beta_map', parameters.surfix,'_abs.nii'];
    spm_write_vol(variables.vo, abs(output.beta_map));	
    
    %% compute the sensitivity map
    alpha = m.sv_coef.';
    SVs = m.SVs;
    X = SVs.';
%     gamma = 1;
    K = zeros(m.totalSV, m.totalSV);
    for ni=1:m.totalSV
        for nj=ni:m.totalSV
            K(ni, nj) = exp(variables.E(m.sv_indices(ni), m.sv_indices(nj)) / parameters.gamma);
            K(nj, ni) = K(ni, nj);
        end
    end
        
    map = X*diag(alpha)*K - X*diag(alpha*K);
    sensitivity_map = sum(map, 2)/numel(alpha);

    img_scale = 10/max(abs(sensitivity_map(:)));
    sensitivity_map = img_scale * sensitivity_map;
    output.sensitivity_map = zeros(nx, ny, nz);
    output.sensitivity_map(variables.l_idx) = sensitivity_map;
    
    variables.vo.fname = ['./', parameters.timestamp, '/', ...
        parameters.score_name, '/SVR_LSM/sensitivity_map', parameters.surfix,'.nii'];
    spm_write_vol(variables.vo, output.sensitivity_map);
%     variables.vo.fname = ['./', parameters.timestamp, '/', ...
%         parameters.score_name, '/SVR_LSM/s_map', parameters.surfix,'_abs.nii'];
%     spm_write_vol(variables.vo, abs(output.pre_img));	
    
    variables.pos_idx = find(output.sensitivity_map>0);
    variables.neg_idx = find(output.sensitivity_map<0);
    
    variables.img_scale = img_scale;
    
    % compute the squared pre-image
    sensitivity_map_squared = sum(map.*map, 2)/numel(alpha);
    
    img_scale2 = 10/max(abs(sensitivity_map_squared(:)));
    sensitivity_map_squared = img_scale2 * sensitivity_map_squared;
    output.sensitivity_map_squared = zeros(nx, ny, nz);
    output.sensitivity_map_squared(variables.l_idx) = sensitivity_map_squared;

    variables.vo.fname = ['./', parameters.timestamp, '/', ...
        parameters.score_name, '/SVR_LSM/s_map', parameters.surfix,'_squared.nii'];
    spm_write_vol(variables.vo, abs(output.sensitivity_map_squared));	
   
    variables.img_scale2 = img_scale2;

end
