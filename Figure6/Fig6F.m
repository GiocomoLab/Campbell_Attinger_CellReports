% Makes the following figures:
% Fig 6F
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

%% Figure 6F right: scatter plot percentage per session

hfig = figure('Position',[200 200 100 400]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6F right: percentage distance tuned neurons scatterplot';
umouse = unique(dt.Mouse);
usession = unique(dt.Session);
plot_col = cbrewer('qual','Set1',numel(umouse));
y = nan(numel(usession),1);
for sIdx = 1:numel(usession)
    in_sesh = strcmp(dt.Session,usession{sIdx});
    mIdx = strcmp(umouse,unique(dt.Mouse(in_sesh)));
    dist_tuned = in_sesh & ...
        dt.pval<opt.pval_cutoff & ...
        dt.prom>opt.min_prom;
    y(sIdx) = sum(dist_tuned)/sum(in_sesh) * 100;
    myscatter(randn(1,1),y(sIdx),plot_col(mIdx,:),0.4);
end
ylabel('% distance tuned neurons');
xlim(xlim*1.5);
h = gca; h.XAxis.Visible = 'off';

%% Figure 6F left: scatter plot num per session

hfig = figure('Position',[200 200 100 400]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6F left: num distance tuned neurons scatterplot';
umouse = unique(dt.Mouse);
usession = unique(dt.Session);
plot_col = cbrewer('qual','Set1',numel(umouse));
y = nan(numel(usession),1);
for sIdx = 1:numel(usession)
    in_sesh = strcmp(dt.Session,usession{sIdx});
    mIdx = strcmp(umouse,unique(dt.Mouse(in_sesh)));
    dist_tuned = in_sesh & ...
        dt.pval<opt.pval_cutoff & ...
        dt.prom>opt.min_prom;
    y(sIdx) = sum(dist_tuned);
    myscatter(randn(1,1),y(sIdx),plot_col(mIdx,:),0.4);
end
ylabel('Num. distance tuned neurons');
xlim(xlim*1.5);
h = gca; h.XAxis.Visible = 'off';
