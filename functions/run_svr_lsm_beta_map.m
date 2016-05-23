function [variables] = run_svr_lsm_beta_map(parameters, variables)
    %% get the beta-map
    
    variables.output_folder = ['./', parameters.timestamp, '/', parameters.score_name, '/SVR_LSM_beta_map', parameters.surfix];
    
    if(~exist(variables.output_folder))
        eval(['mkdir ', variables.output_folder]);
    end
    
    cmd = ['-s 3 -t 2 -c ', num2str(parameters.cost), ' -g ', num2str(parameters.gamma), ' -q'];
    variables.one_score = variables.one_score*100/max(abs(variables.one_score));

    [beta_map, variables] = get_beta_map(parameters, variables, cmd);

    variables.atlas_aal_roi_filename = [variables.output_folder, '/beta_map_aal.txt'];
    variables.atlas_ba_roi_filename = [variables.output_folder, '/beta_map_brodmann.txt'];
    [atlas_aal_roi, atlas_ba_roi] = atlas_vol_analysis(parameters, variables, ...
        beta_map, quantile(abs(beta_map(variables.m_idx)), parameters.percentile/100));

    %% Permutation test
    if(parameters.PermNum)
        variables.output_folder = [variables.output_folder, '/PMU_Iter_', num2str(parameters.PermNum)];
        if(~exist([variables.output_folder]))
            eval(['mkdir ', variables.output_folder]);
        end
    
        [p_map, p_map_pos, p_map_neg, variables] = run_beta_PMU(parameters, variables, cmd, beta_map);

        if(variables.q_fdr>0)
            variables.atlas_aal_roi_filename = [variables.output_folder, '/PMU_map_aal.txt'];
            variables.atlas_ba_roi_filename = [variables.output_folder, '/PMU_map_brodmann.txt'];
            [atlas_aal_roi, atlas_ba_roi] = PMU_atlas_vol_analysis(parameters, variables, p_map);
        end
    end
end