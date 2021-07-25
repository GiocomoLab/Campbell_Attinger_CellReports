% Makes the following figures:
% Fig 6D
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = fullfile(paths.intermediate_data,'dark_distance');

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% load dist tuning
binsize = 2;
smoothsig = 4;
dt = load(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',binsize,smoothsig)),'dist_tuning','opt');
opt = dt.opt;
dt = dt.dist_tuning;

% analysis options
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence

dist_tuned = dt.pval<opt.pval_cutoff & dt.prom>opt.min_prom;

%% compute peak zscore
peak_zscore = nan(size(dt,1),1);
for cIdx = 1:size(dt,1)
    peak_zscore(cIdx) = (dt.peak(cIdx)-nanmean(dt.peak_shuf(cIdx,:)))/nanstd(dt.peak_shuf(cIdx,:));
end

%% Figure 6D: plot prom vs peak zscore
hfig = figure('Position',[300 300 350 300]); hold on;
hfig.Name = 'Figure 6D: prom vs peak zscore';
myscatter(dt.prom(dist_tuned),peak_zscore(dist_tuned),lines(1),0.2);
myscatter(dt.prom(~dist_tuned),peak_zscore(~dist_tuned),'k',0.2);
xlabel('Peak prominence');
ylabel('Peak z-score vs. shuffle');
legend({sprintf('Dist (%d)',sum(dist_tuned)),sprintf('Non-Dist (%d)',sum(~dist_tuned))},'Location','southeast');
set_fig_prefs;
xlim([-0.02 0.8]);
yticks(-10:10:30);