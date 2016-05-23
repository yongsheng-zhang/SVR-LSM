function [roi_aal, roi_ba] = get_significant_roi(variables, mask_idx)
    if(length(mask_idx) == 0)
        roi_aal = [];
        roi_ba = [];
        return;
    end
    %% aal atlas
    ba_value = variables.atlas_aal_map(mask_idx);
    ba_min = min(ba_value);
    ba_max = max(ba_value);
    
    count = 0;
    ba_idx = 0;
    ba_count = 0;
    for ni = 1:ba_max
        tmp_idx = find(ba_value == ni);
        ba_idx(ni,1) = ni;
        count(ni,1) = length(tmp_idx);
        ba_count(ni,1) = length(find(variables.atlas_aal_map(:)==ni));
    end

    na_idx = find(count == 0);
    ba_idx(na_idx) = [];
    count(na_idx) = [];
    ba_count(na_idx) = [];

    if(length(ba_idx) == 0)
        roi_aal = [];
        roi_ba = [];
        return;
    end
    roi_aal(:,1) = double(ba_idx);
    roi_aal(:,2) = double(count);
    roi_aal(:,3) = double(ba_count);
    
    %% Brodmann's Area atlas
    ba_value = variables.atlas_ba_map(mask_idx);
    ba_min = min(ba_value);
    ba_max = max(ba_value);

    count = 0;
    ba_idx = 0;
    ba_count = 0;
    for ni = 1:ba_max
        tmp_idx = find(ba_value == ni);
        ba_idx(ni,1) = ni;
        count(ni,1) = length(tmp_idx);
        ba_count(ni,1) = length(find(variables.atlas_ba_map(:)==ni));
    end

    na_idx = find(count == 0);
    ba_idx(na_idx) = [];
    count(na_idx) = [];
    ba_count(na_idx) = [];
    
    if(length(ba_idx)==0)
        roi_ba(:,1) = 0;
        roi_ba(:,2) = 0;
        roi_ba(:,3) = 0;
    else
        roi_ba(:,1) = double(ba_idx);
        roi_ba(:,2) = double(count);
        roi_ba(:,3) = double(ba_count);
    end
end