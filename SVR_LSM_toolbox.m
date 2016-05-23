% SVR-based multivariate lesion-symptom mapping script
% Created by Yongsheng Zhang, 2014
% yongsheng.zhang.cn@gmail.com
% Please cite the following paper if you think this script is helpful:
%
% Yongsheng Zhang, Daniel Y. Kimberg, H. Branch Coslett, Myrna F. Schwartz,
% Ze Wang, "Multivariate Lesion-Symptom Mapping Using Support Vector
% Regression", Human Brain Mapping, 2014.

clear
%% add path of SPM, libsvm and support functions of SVR-LSM
addpath('/pkg/spm12');% add spm folder to search path
addpath('./libsvm-3.18/matlab'); % add libSVpuM's runtime folder to search path
addpath('./functions'); % folder for support functions

%% setup batch parameters
parameters.score_file = './score/PNT.csv';
parameters.score_name = 'Sim_ROI_123';
parameters.control_variable_name = ''; % control variable filename
parameters.lesion_img_folder = './lesion_imgs/'; % file or directory that contain lesion image, 

% SVR-LSM parameters
parameters.dtlvc_flag = 1; % 1 for dTLVC for SVR-LSM; 0 without dTLVC
parameters.svr_map_type = 1; % 1 for beta-map, 2 for sensitivity map, 3 for both
parameters.invert_p_map_flag = 1; % Inverse p-map, i.e., use 1-p instead of p for display on MRIcron.
parameters.cost = 30;
parameters.gamma = 2;

% permutation test parameters
parameters.PermNum = 10; % iteration times for permutation test, 0 for no permutation test
parameters.lesion_thresh = 10; % The least lesion subject number for a voxel to be considered in the following analysis. 

% atlas file name for calculate voxel numbers in each atlas region.
parameters.q_thresh = 0.05;
parameters.percentile = 95; 
parameters.min_cluster_size = 50;
parameters.atlas_aal = './atlas/aal.nii';
parameters.atlas_aal_name = './atlas/aal.nii.txt'; % store the brain area name for A atlas
parameters.atlas_brodmann = './atlas/brodmann.nii';

%% print authorship information
clc
fprintf('\n%s',repmat(sprintf('='),1,40))
begin_time = clock;
fprintf(' %d-%d-%d %d:%d:%.0f\n',begin_time(1), begin_time(2), begin_time(3), ...
    begin_time(4), begin_time(5), begin_time(6));
fprintf('Running multivariate lesion-symptom mapping.\n')
fprintf('More details about the method, please refer to :\n');
fprintf('Yongsheng Zhang, Daniel Y. Kimberg, H. Branch Coslett, Myrna F. Schwartz,  Ze Wang, \n');
fprintf('\"Multivariate Lesion-Symptom Mapping Using Support Vector Regression\", \n');
fprintf('Human Brain Mapping, 2014.\n');
fprintf('\n%s\n',repmat(sprintf('-'),1,40))

parameters.timestamp = date;
variables.output_folder = ['./', parameters.timestamp, '/', parameters.score_name];

if(~exist(variables.output_folder))
    eval(['mkdir ', variables.output_folder]);
end

copyfile([mfilename, '.m'], fullfile(variables.output_folder, 'MyConfig.bak'));

if(and(length(parameters.control_variable_name), parameters.dtlvc_flag))
    error('Only one total lesion volume control method should be used, either set parameters.dtlvc_flag to 0 or set parameters.control_variable_name to a blank string!');
end
%% Read behavior score and lesion images
[variables] = read_behavior_score(parameters);
[variables] = read_lesion_imgs(parameters, variables);

%% regression out control variables from the behavior score of interest
fprintf('Running analysis for ''%s'' ', parameters.score_name);
if(~isempty(parameters.control_variable_name))
    parameters.surfix = ['_Ctrl_', parameters.control_variable_name];
    variables.one_score = regressout_variables(variables.one_score, variables.control_var);
    fprintf('by regressing out variable: ''%s''. \n', parameters.control_variable_name);
else
    parameters.surfix = '';
    fprintf('without variable control. \n');
end

%% preprocessing lesion-map if dTLVC is desired
if parameters.dtlvc_flag
    parameters.surfix = [parameters.surfix, '_dTLVC'];    
    fprintf('lesion volume: \n');
    for ni=1:variables.SubNum
        variables.lesion_dat(ni,:) = variables.lesion_dat(ni,:)/sqrt(variables.lesion_vol(ni,1));
    end
end

%% SVR-LSM beta-map
if(rem(parameters.svr_map_type,2))
    [variables] = run_svr_lsm_beta_map(parameters, variables);
end

%% SVR-LSM sensitivity-map
if(uint8(floor(parameters.svr_map_type/2)))
    fprintf('\n\nPlease refer to the following paper for details about sensitivity map:\n');
    fprintf('Rasmussen PM, Madsen KH, Lund TE, Hansen LK (2011): \n');
    fprintf('Visualization of nonlinear kernel models in neuroimaging by sensitivity maps.\n');
    fprintf('NeuroImages 55:1120-1131\n');
    [variables] = run_svr_lsm_sensitivity_map(parameters, variables);
end
