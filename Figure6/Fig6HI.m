% Makes the following figures:
% Fig 6H,I
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

%% Figure 6H: all together

hfig = figure('Position',[200 200 400 200]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6H: preferred distance histogram all together';
y = dt.period(dist_tuned);
distbinedge = 0:2:opt.max_lag;
distbincent = distbinedge(1:end-1)+mean(diff(distbinedge))/2;
xid = distbincent;
hcounts = histcounts(y,distbinedge);
%hcounts = interp1(distbincent,hcounts,xid);
hcounts = gauss_smoothing(hcounts,1);
plot(xid,hcounts);
[~,x] = findpeaks(hcounts);
x = xid(x(1));
x = 75;
for i = 1:7
    plot(repmat(x*sqrt(2)^(i-1),1,2),ylim,'k:');
end
xlabel('Preferred distance (cm)')
ylabel('Num. dist. cells');
title('All mice');

%% Figure 6I: by mouse

hfig = figure('Position',[200 200 400 200]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6I: preferred distance histogram by mouse';
umouse = unique(dt.Mouse);
plot_col = cbrewer('qual','Set1',numel(umouse));
for mIdx = 1:numel(umouse)
    keep = strcmp(dt.Mouse,umouse{mIdx}) & dist_tuned;
    y = dt.period(keep);
    hcounts = histcounts(y,distbinedge);
    %hcounts = interp1(distbincent,hcounts,xid);
    hcounts = gauss_smoothing(hcounts,1);
    plot(xid,hcounts,'Color',plot_col(mIdx,:));
end
xlabel('Preferred distance (cm)');
ylabel('Num. dist. cells');
title('By mouse');
leg_mouse = umouse;
for mIdx = 1:numel(leg_mouse)
    leg_mouse{mIdx} = leg_mouse{mIdx}(3:4);
end
legend(leg_mouse);
for i = 1:7
    plot(repmat(x*sqrt(2)^(i-1),1,2),ylim,'k:','HandleVisibility','off');
end