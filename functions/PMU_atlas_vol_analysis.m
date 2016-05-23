function [atlas_aal_roi, atlas_ba_roi] = PMU_atlas_vol_analysis(parameters, variables, p_map)
    inv_q_map = zeros(size(p_map));
    inv_q_map(variables.l_idx) = 1-p_map(variables.l_idx);
    q_idx = find(inv_q_map <= 1-variables.q_fdr);
    inv_q_map(q_idx) = 0;
    inv_q_map = remove_scatter_clusters(inv_q_map, parameters.min_cluster_size);
    mask_idx = find(inv_q_map > 1-variables.q_fdr);
%     inv_p_map = zeros(size(p_map));
%     inv_p_map(variables.l_idx) = 1-p_map(variables.l_idx);
%     mask_idx = find(inv_p_map > 1-variables.q_fdr); 
    [atlas_aal_roi, atlas_ba_roi] = get_significant_roi(variables, mask_idx);

    % output the AAL atlas area name and significant voxel number for highlited ROIs
    fid = fopen(variables.atlas_aal_roi_filename, 'w');
    fprintf(fid, '============ FDR q < %7.4f, equivalent to p < %7.4f \n', parameters.q_thresh, variables.q_fdr);
    if(size(atlas_aal_roi,1) == 0)
        fprintf(fid, '\nNone\n');
    else
        fprintf(fid, 'index | %15s    | \t voxel number\n', 'atlas_name');
        for ni=1:size(atlas_aal_roi,1)
            fprintf(fid, ' %3d, %20s, %6d out of %6d voxels\n', atlas_aal_roi(ni,1), variables.atlas_aal.name{atlas_aal_roi(ni,1)}, atlas_aal_roi(ni,2), atlas_aal_roi(ni,3));
        end
    end
    if(~isempty(variables.exclude_idx))
        fprintf(fid, '*******************************************************************\n');
        fprintf(fid, 'Warning: The following subjects have been excluded from analysis because they have no survival voxel after thresholding, \n');
        for(ni=1:length(variables.exclude_idx))
            fprintf(fid, '%s \n', variables.excluded_SubjectID{ni});
        end
        
    end    
    fclose(fid);

    % output the Brodmann Area name and significant voxel number for highlited ROIs
    fid = fopen(variables.atlas_ba_roi_filename, 'w');
    fprintf(fid, '============ FDR q < %7.4f, equivalent to p < %7.4f \n', parameters.q_thresh, variables.q_fdr);
    if(size(atlas_ba_roi,1) == 0)
        fprintf(fid, '\nNone\n');
    else
        fprintf(fid, 'index | voxel number\n', 'atlas_name');
        for ni=1:size(atlas_ba_roi,1)
            fprintf(fid, ' %3d, %6d out of %6d voxels\n', atlas_ba_roi(ni,1), atlas_ba_roi(ni,2), atlas_ba_roi(ni,3));
        end
    end
    if(~isempty(variables.exclude_idx))
        fprintf(fid, '*******************************************************************\n');
        fprintf(fid, 'Warning: The following subjects have been excluded from analysis because they have no survival voxel after thresholding, \n');
        for(ni=1:length(variables.exclude_idx))
            fprintf(fid, '%s \n', variables.excluded_SubjectID{ni});
        end
        
    end    
    fclose(fid);        
end
