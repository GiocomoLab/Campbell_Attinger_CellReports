% preps workspace paths and variables
% MGC 10/1/2019

%% paths
data_dir = fullfile(neuropix_folder,'data');
cell_info_dir = fullfile(neuropix_folder,'matlab_scripts','cell_info');
session_info_dir = fullfile(neuropix_folder,'matlab_scripts','session_info');
cell_info_file = 'cell_info_ALL_Oct2019_v2';

% make folders if they don't exist already
data_save_dir = fullfile(neuropix_folder,data_save_dir);
image_save_dir = fullfile(neuropix_folder,image_save_dir);
if exist('data_save_dir','var')
    if exist(data_save_dir,'dir')~=7
        mkdir(data_save_dir);
    end
end
if exist('image_save_dir','var')
    if exist(image_save_dir,'dir')~=7
        mkdir(image_save_dir);
    end
end

%% params that mostly stay fixed
tol = 0.0001; % for comparing floats
num_trials_per_stability_block = 4;
params = readtable(fullfile(neuropix_folder,'UniversalParams.xlsx'));
numtrialspergainchange = 4;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;
track_length = params.TrackEnd-params.TrackStart;

%% load cell_info table
load(fullfile(cell_info_dir,cell_info_file));

%% load session_file
if exist('session_file','var')
    load(fullfile(session_info_dir,session_file));
    mouse = {};
    mouse_date = {};
    for i = 1:numel(session_name)
        mouse_date_this = strsplit(session_name{i},'_');
        mouse{i} = mouse_date_this{1};
        mouse_date{i} = strcat(mouse_date_this{1},'_',mouse_date_this{2});
    end
    mouse = mouse';
    mouse_uniq = unique(mouse);
    mouse_date = mouse_date';
    mouse_date_uniq = unique(mouse_date);
    clear mouse_date_this i
end