function [beta_map, variables] = get_beta_map(parameters, variables, cmd)

    fprintf('SVR-LSM beta map...');
    
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
    %% compute the beta-map
    w = m.sv_coef.'*m.SVs;
    variables.beta_scale = 10/max(abs(w));
    tmp = zeros(nx, ny, nz);
    tmp(variables.l_idx) = w.'*variables.beta_scale;
    beta_map = zeros(nx, ny, nz);
    beta_map(variables.m_idx) = tmp(variables.m_idx);
    
    variables.vo.fname = [variables.output_folder, '/beta_map.nii'];
    spm_write_vol(variables.vo, beta_map);
    variables.vo.fname = [variables.output_folder, '/beta_map_abs.nii'];
    spm_write_vol(variables.vo, abs(beta_map));	
    
    top_beta_map = beta_map;
    non_sig_idx = find(abs(top_beta_map) <= quantile(abs(beta_map(variables.m_idx)), ...
        parameters.percentile/100));    
    top_beta_map(non_sig_idx) = 0;
    top_beta_map = remove_scatter_clusters(top_beta_map, parameters.min_cluster_size);
    variables.vo.fname = [variables.output_folder, '/beta_map_threshold_', ...
        num2str(quantile(abs(beta_map(variables.m_idx)), parameters.percentile/100)),'.nii'];
    spm_write_vol(variables.vo, abs(top_beta_map));
    
    variables.pos_idx = find(beta_map>0);
    variables.neg_idx = find(beta_map<0);

end
