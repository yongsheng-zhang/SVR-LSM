function [sensitivity_map, variables] = get_sensitivity_map(parameters, variables, cmd)

    fprintf('SVR-LSM sensitivity map...');
    
    %%
    SVs = variables.lesion_dat.';

    E = zeros(variables.SubNum, variables.SubNum);
    for ni=1:variables.SubNum
        for nj=ni:variables.SubNum
            E(ni, nj) = -1*norm(SVs(ni,:) - SVs(nj,:))^2;
            E(nj, ni) = E(ni, nj);
        end
    end
    variables.E = E;
    %%
    nx = variables.vo.dim(1);
    ny = variables.vo.dim(2);
    nz = variables.vo.dim(3);
    
    if(size(variables.one_score,1)<size(variables.one_score,2))
        variables.one_score = variables.one_score.';
    end
    m = svmtrain(variables.one_score,sparse(variables.lesion_dat),cmd);
    
    fprintf('done.\n');
    fprintf('# of support vectors: %d\n', m.totalSV);
    fprintf('# of bounded support vectors: %d\n', sum(m.sv_coef == parameters.cost));
    %% compute the sensitivity map
    alpha = m.sv_coef.';
    SVs = m.SVs;
    X = SVs.';
    K = zeros(m.totalSV, m.totalSV);
    for ni=1:m.totalSV
        for nj=ni:m.totalSV
            K(ni, nj) = exp(variables.E(m.sv_indices(ni), m.sv_indices(nj)) / parameters.gamma);
            K(nj, ni) = K(ni, nj);
        end
    end
        
    map = X*diag(alpha)*K - X*diag(alpha*K);
    s_map = sum(map, 2)/numel(alpha);

    variables.s_scale = 10/max(abs(s_map(:)));
    sensitivity_map = variables.s_scale * s_map;
    sensitivity_map = zeros(nx, ny, nz);
    sensitivity_map(variables.l_idx) = s_map*variables.s_scale;
    
    variables.vo.fname = [variables.output_folder, '/sensitivity_map.nii'];
    spm_write_vol(variables.vo, sensitivity_map);
    variables.vo.fname = [variables.output_folder, '/sensitivity_map_abs.nii'];
    spm_write_vol(variables.vo, abs(sensitivity_map));	
    
    top_s_map = sensitivity_map;
    non_sig_idx = find(abs(top_s_map) <= quantile(abs(sensitivity_map(variables.l_idx)), ...
        parameters.percentile/100));    
    top_s_map(non_sig_idx) = 0;
    top_s_map = remove_scatter_clusters(top_s_map, parameters.min_cluster_size);
    variables.vo.fname = [variables.output_folder, '/sensitivity_map_threshold_', ...
        num2str(quantile(abs(sensitivity_map(variables.l_idx)), parameters.percentile/100)),'.nii'];
    spm_write_vol(variables.vo, abs(top_s_map));    
    
    variables.pos_idx = find(sensitivity_map>0);
    variables.neg_idx = find(sensitivity_map<0);   

end
