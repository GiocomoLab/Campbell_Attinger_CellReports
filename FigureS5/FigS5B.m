% Makes the following figures:
% Fig S5B
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

%% Figure S5B, left: by mouse, subplot

leg_mouse = unique(dt.Mouse);
for mIdx = 1:numel(leg_mouse)
    leg_mouse{mIdx} = leg_mouse{mIdx}(3:4);
end

hfig = figure('Position',[200 200 1200 800]); set_fig_prefs;
hfig.Name = 'Figure S5B, left: preferred distance histogram by mouse subplot';
umouse = unique(dt.Mouse);
distbinedge_this = 0:4:800;
xid_this = 2:4:798;
ratio = [];
for mIdx = 1:numel(umouse)
    subplot(4,3,mIdx);
    keep = strcmp(dt.Mouse,umouse{mIdx}) & dist_tuned;
    y = dt.period(keep);
    hcounts = histcounts(y,distbinedge_this);
    %hcounts = interp1(distbincent,hcounts,xid);
    hcounts = gauss_smoothing(hcounts,1);
    [~,locs,~,prom] = findpeaks(hcounts);
    keep_pks = prom>max(hcounts)/8;
    locs = locs(keep_pks);
    if numel(locs)>1
        ratio = [ratio locs(2:end)./locs(1:end-1)];
    end
    stairs(xid_this,hcounts); hold on;
    plot(xid_this(locs),hcounts(locs),'ro','MarkerFaceColor','r');
    title(leg_mouse{mIdx});
    box off;
    xlabel('Preferred distance (cm)');
end


%% Figure S5B right: histograms of module ratios

hfig = figure('Position',[400 400 400 300]); hold on;
hfig.Name = 'Figure S5B, right: histogram of module ratios';
hbar=histogram(ratio,1:0.05:2.5,'EdgeColor','k','FaceColor','k','FaceAlpha',1);
hbar.HandleVisibility='off';
box off;
yticks(0:5);
plot([sqrt(2) sqrt(2)],[0 5],'g');
plot([mean(ratio) mean(ratio)],[0 5],'r');
legend('sqrt(2)','mean');
ylabel('Count');
xlabel('Module ratio');