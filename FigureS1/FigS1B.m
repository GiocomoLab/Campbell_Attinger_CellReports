% Makes the following figures:
% Fig S1B left
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

%restoredefaultpath;

% add helper functions to path
%addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***

% *** NOTE: Download spikes repository from https://github.com/cortex-lab/spikes
addpath(genpath(paths.spikes_repo));

% load cell_info table
load(fullfile(paths.intermediate_data,fullfile('cell_info','cell_info_Apr2020'))); 

opt = load_default_opt;
opt.plot_col = cbrewer('qual','Set2',3);

session_this = 'npJ3_0506_contrast_1';

%%
dat = load(fullfile(paths.data,session_this));

good_cells = dat.sp.cids(dat.sp.cgs==2);
cellid_in_br = cell_info.CellID(strcmp(cell_info.Session,session_this) & strcmp(cell_info.BrainRegion,'MEC'));

in_br = ismember(dat.sp.clu,cellid_in_br);

[~,depth_all] = templatePositionsAmplitudes(dat.sp.temps,dat.sp.winv,dat.sp.ycoords,dat.sp.spikeTemplates,dat.sp.tempScalingAmps);
%depth_all = dat.anatomy.tip_distance;
max_depth_all = max(depth_all(ismember(dat.sp.clu,good_cells)));
min_depth_all = min(depth_all(ismember(dat.sp.clu,good_cells)));
max_depth_in_br = max(depth_all(in_br));
min_depth_in_br = min(depth_all(in_br));


%% Figure S1B left
hfig = figure('Position',[300 300 800 400]); hold on;
hfig.Name = sprintf('Figure S1B left: MEC %s t_start=1327 with inset',session_this);

start_time = 1327;
stop_time = 1331;

keep_spike = dat.sp.st>start_time & dat.sp.st<stop_time & depth_all<=max_depth_all & depth_all>=min_depth_all;
in_br2 = depth_all>=min_depth_in_br & depth_all<=max_depth_in_br;
in_br2 = in_br2(keep_spike);

spiket_snippet = dat.sp.st(keep_spike)-start_time;
depth_snippet = depth_all(keep_spike);


subplot(1,3,1:2); hold on
plot(spiket_snippet(in_br2),depth_snippet(in_br2),'.','Color',opt.plot_col(1,:));
plot(spiket_snippet(~in_br2),depth_snippet(~in_br2),'k.');
xlim([0 4]);
ylim([min(depth_snippet)-100 max(depth_snippet)+100]);
xlabel('time (sec)');
ylabel('depth from probe tip');
title(sprintf('%s (MEC)',session_this),'Interpreter','none');
ymin = min(ylim);
patch([2 2.2 2.2 2],[ymin ymin 1700 1700],'r','EdgeColor','r','FaceColor','none'); 
box off;


% inset

start_time = 1329;
stop_time = 1329.2;

keep_spike = dat.sp.st>start_time & dat.sp.st<stop_time & depth_all<=max_depth_all & depth_all>=min_depth_all;
in_br2 = depth_all>=min_depth_in_br & depth_all<=max_depth_in_br;
in_br2 = in_br2(keep_spike);

spiket_snippet = dat.sp.st(keep_spike)-start_time;
depth_snippet = depth_all(keep_spike);

subplot(1,3,3); hold on;
plot(spiket_snippet(in_br2),depth_snippet(in_br2),'.','Color',opt.plot_col(1,:));
plot(spiket_snippet(~in_br2),depth_snippet(~in_br2),'k.');
xlim([0 .2]);
ylim([min(depth_snippet)-100 1700]);
xlabel('time (sec)');
box off;
yticks([]);