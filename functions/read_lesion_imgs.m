function [variables] = read_lesion_imgs(parameters, variables)

%% Read the lesion data and get the lesion volume of each subject
variables.lesion_vol = zeros(size(variables.one_score));

for ni=1:length(variables.SubjectID)
    fname = fullfile(parameters.lesion_img_folder, [variables.SubjectID{ni}, '.nii']);
    if(~exist(fname))
        warning('can not find lesion image file: %s\n', fname);
        continue;
    end
    vo = spm_vol(fname);
    tmp = spm_read_vols(vo);
    
    tmp(find(isnan(tmp(:))==1)) = 0;
    Ldat(:,:,:,ni) = uint8(tmp);
    
    variables.lesion_vol(ni,1) = sum(tmp(:));
end
variables.vo = vo;
variables.vo.name = 'NULL.nii';
variables.vo.dt = [64,0];

%% get a mask based on overlapping map and given mask image
mask_map = sum(Ldat, 4);
% index of voxels with lesion on at least 1 subject
variables.l_idx = find(mask_map >= 1); 
% index of voxels that will be excluded from result
variables.m_idx = find(mask_map >= parameters.lesion_thresh); 

% get the voxel index within the generated mask
Mdat = reshape(Ldat, length(mask_map(:)), variables.SubNum).';

variables.lesion_dat = double(Mdat(:,variables.l_idx));

% Check the thresholded lesion data, remove subjects who have no
% survival voxel within the mask.
Ldat_tmp = double(Mdat(:,variables.m_idx));
sub_idx = find(sum(Ldat_tmp, 2) == 0);
variables.exclude_idx = sub_idx;
variables.SubNum = variables.SubNum - length(sub_idx);
if(~isempty(sub_idx))
    fprintf('<strong> Warning: The following subjects have no survival voxel after thresholding, \nand so will be excluded in the following analysis:\n</strong>')
    for(ni=1:length(sub_idx))
        fprintf('<strong>%s </strong>\n', variables.SubjectID{sub_idx(ni)})
        variables.lesion_dat(ni,:) = [];
        variables.one_score(ni) = [];
        variables.lesion_vol(ni) = [];
        variables.excluded_SubjectID{ni} = variables.SubjectID{sub_idx(ni)};
        variables.SubjectID(sub_idx(ni)) = [];
        
    end
    fprintf('\n')
end
%% Read atlas images
vo = spm_vol(parameters.atlas_aal);
variables.atlas_aal_map = uint8(spm_read_vols(vo));
[variables.atlas_aal.id, variables.atlas_aal.name, ba_tmp] = textread(parameters.atlas_aal_name, '%d %s %d');

vo = spm_vol(parameters.atlas_brodmann);
variables.atlas_ba_map = uint8(spm_read_vols(vo));

%% initialize some parameters
variables.neg_idx = [];
variables.pos_idx = [];

end