% Makes the following figures:
% Fig 6G
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = paths.intermediate_data;

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% load dist tuning
paths.rez = fullfile(paths.rez,'dark_distance');
binsize = 2;
smoothsig = 4;
dt = load(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',binsize,smoothsig)),'dist_tuning','opt');
opt = dt.opt;
dt = dt.dist_tuning;

% analysis options
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence

dist_tuned = dt.pval<opt.pval_cutoff & dt.prom>opt.min_prom;

%% Figure 6G: scatter plot of dist neuron period vs depth

x = dt.period(dist_tuned);
y = dt.DepthAdjusted(dist_tuned);
y_all = dt.DepthAdjusted;
xedges = 0:5:800;
yedges = -2000:50:1000;
xcent = xedges(1:end-1)+mean(diff(xedges))/2;
ycent = yedges(1:end-1)+mean(diff(yedges))/2;

hfig = figure('Position',[300 300 500 300]);
hfig.Name = 'Figure 6G: scatter of dist neuron period vs depth';

subplot(1,4,1:3); hold on;
myscatter(x,y,lines(1),0.15,15);
xticks(0:200:800);
yticks(-2000:1000:1000);
yticklabels(-2:1);
xlabel('Period (cm)');
ylabel('Dist. from MEC entry (mm)');
set(gca,'TickDir','out');
box off;
ylim([-2000 1000]);
plot(xlim,[0 0],'k--');

subplot(1,4,4); hold on; set_fig_prefs;
n1 = histcounts(y,yedges);
n2 = histcounts(y_all,yedges);
plot(100*n1./n2,ycent,'Color',lines(1));
box off
plot(xlim,[0 0],'k--');
yticklabels([]);
xlabel('% Dist. Cells');