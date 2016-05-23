function [p_map, p_map_pos, p_map_neg, variables] = run_beta_PMU(parameters, variables, cmd, beta_map)
    
    fprintf('\nBegin to run SVR-LSM Permutation Test...\n');

    nx = variables.vo.dim(1);
    ny = variables.vo.dim(2);
    nz = variables.vo.dim(3);
    
    map_count = zeros(1,length(variables.m_idx)); 
    map_count_pos = zeros(1,length(variables.m_idx));
    map_count_neg = zeros(1,length(variables.m_idx));
    beta_pos_idx = find(beta_map(variables.m_idx)>=0);
    beta_neg_idx = find(beta_map(variables.m_idx)<0);
    ori_beta_val = beta_map(variables.m_idx).';
        
    tic
    prompt_str = '';
    for PermIdx=1:parameters.PermNum
        % random permute subjects order
        loc = randperm(length(variables.one_score));
        trial_score = variables.one_score(loc);

        m = svmtrain(trial_score,sparse(variables.lesion_dat),cmd);
        % compute the beta map
        alpha = m.sv_coef.';
        SVs = m.SVs;
        pmu_beta_map = variables.beta_scale * alpha*SVs;  
        
        tmp_map = zeros(nx, ny, nz);
        tmp_map(variables.l_idx) = pmu_beta_map;
        pmu_beta_map = tmp_map(variables.m_idx).';

        map_count = map_count + (abs(pmu_beta_map) > abs(ori_beta_val));
        map_count_pos(beta_pos_idx) = map_count_pos(beta_pos_idx) + (pmu_beta_map(beta_pos_idx) > ori_beta_val(beta_pos_idx));
        map_count_neg(beta_neg_idx) = map_count_neg(beta_neg_idx) + (pmu_beta_map(beta_neg_idx) < ori_beta_val(beta_neg_idx));

        elapsed_time = toc;
        remain_time = round(elapsed_time * (parameters.PermNum - PermIdx)/(PermIdx));
        remain_time_h = floor(remain_time/3600);
        remain_time_m = floor((remain_time - remain_time_h*3600)/60);
        remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);

        if(~isempty(prompt_str))
            fprintf(repmat('\b', size(prompt_str)));
        end
        prompt_str = sprintf(['Permutation-', num2str(PermIdx), '/', num2str(parameters.PermNum), ...
            ': ', datestr(clock), '\tEstimated remaining time: ', num2str(remain_time_h), ' h ', ...
            num2str(remain_time_m), ' m ' num2str(remain_time_s), 's\n']);

        fprintf(prompt_str);
    end
    
    p_map = zeros(nx, ny, nz);
    p_map(variables.m_idx) = map_count/(parameters.PermNum+1);
    p_map_pos = zeros(nx, ny, nz);
    p_map_pos(variables.m_idx) = map_count_pos/(parameters.PermNum+1);
    p_map_neg = zeros(nx, ny, nz);
    p_map_neg(variables.m_idx) = map_count_neg/(parameters.PermNum+1);     

    if(parameters.invert_p_map_flag)
        % PMU result for beta-map
        tmp_out_img = zeros(nx, ny, nz);
        tmp_out_img(variables.m_idx) = 1-p_map(variables.m_idx);
        variables.vo.fname = [variables.output_folder, '/p_map_beta_inverse.nii'];
        spm_write_vol(variables.vo, tmp_out_img);

        tmp_out_img = zeros(nx, ny, nz);
        tmp_out_img(variables.pos_idx) = 1-p_map_pos(variables.pos_idx);
        variables.vo.fname = [variables.output_folder, '/p_map_beta_pos_inverse.nii'];
        spm_write_vol(variables.vo, tmp_out_img);                   

        tmp_out_img = zeros(nx, ny, nz);
        tmp_out_img(variables.neg_idx) = 1-p_map_neg(variables.neg_idx);                
        variables.vo.fname = [variables.output_folder, '/p_map_beta_neg_inverse.nii'];
        spm_write_vol(variables.vo, tmp_out_img);   

    else
        variables.vo.fname = [variables.output_folder, '/p_map.nii'];
        spm_write_vol(variables.vo, p_map);        
        variables.vo.fname = [variables.output_folder, '/p_map_pos.nii'];
        spm_write_vol(variables.vo, p_map_pos);            
        variables.vo.fname = [variables.output_folder, '/p_map_neg.nii'];
        spm_write_vol(variables.vo, p_map_neg);    

    end

   
    %% FDR correction
    p_val_roi = map_count/(parameters.PermNum+1);
	[q_fdr, q_masked] = fdr(p_val_roi, parameters.q_thresh);
    variables.q_fdr = q_fdr;
    
    q_map_fdr = zeros(nx, ny, nz);
    if(q_fdr>0)
        q_val_roi = zeros(size(p_val_roi));
        q_idx = find(p_val_roi < q_fdr);
        q_val_roi(q_idx) = 1-p_val_roi(q_idx);
        q_map_fdr(variables.m_idx) = q_val_roi;
        variables.vo.fname = [variables.output_folder, '/q_map_fdr_', ...
            num2str(parameters.q_thresh), '.nii'];
        q_map_fdr = remove_scatter_clusters(q_map_fdr, parameters.min_cluster_size);
        spm_write_vol(variables.vo, q_map_fdr);
    else
        fclose(fopen([variables.output_folder, '/no_voxel_survives_with_FDR_q_', num2str(parameters.q_thresh)], 'w'));
    end

end
