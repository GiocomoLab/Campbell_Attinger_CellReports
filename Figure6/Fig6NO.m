% Makes the following figures:
% Fig 6N,O
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
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1;
opt.min_num_dist_cells = 5; % for by session plot
dt = dt.dist_tuning;
dist_tuned = dt.pval<opt.pval_cutoff & dt.prom>opt.min_prom;

% load sensory delay
load(fullfile(paths.intermediate_data,'sensory_delay','sensory_delay.mat'));

%% compare sensory delay factors between distance cells and non-dist cells

dcells = unique(dt.UniqueID(dist_tuned));
ndcells = unique(dt.UniqueID(~dist_tuned));
assert(isempty(intersect(dcells,ndcells)));

sd_uid = sensory_delay.UniqID;
for cIdx = 1:numel(sd_uid)
    split_this = strsplit(sd_uid{cIdx},'_');
    sd_uid{cIdx} = sprintf('%s_%s_%s',split_this{1},split_this{2},split_this{end});
end

sd_uid_dcells = sd_uid(ismember(sd_uid,dcells));
sd_dcells = sensory_delay.DelayFactor(ismember(sd_uid,dcells));
sd_uid_ndcells = sd_uid(ismember(sd_uid,ndcells));
sd_ndcells = sensory_delay.DelayFactor(ismember(sd_uid,ndcells));

sd_uid_dcells = sd_uid_dcells(~isnan(sd_dcells));
sd_dcells = sd_dcells(~isnan(sd_dcells));
sd_uid_ndcells = sd_uid_ndcells(~isnan(sd_ndcells));
sd_ndcells = sd_ndcells(~isnan(sd_ndcells));

% average by cell
dcells_this = unique(sd_uid_dcells);
ndcells_this = unique(sd_uid_ndcells);
sd_dcells2 = nan(numel(dcells_this),1);
sd_ndcells2 = nan(numel(ndcells_this),1);
for i = 1:numel(dcells_this)
    sd_dcells2(i) = nanmean(sd_dcells(strcmp(sd_uid_dcells,dcells_this{i})));
end
for i = 1:numel(ndcells_this)
    sd_ndcells2(i) = nanmean(sd_ndcells(strcmp(sd_uid_ndcells,ndcells_this{i})));
end
sd_uid_dcells = dcells_this;
sd_uid_ndcells = ndcells_this;

% flipped signs 7 Aug 2020
sd_dcells = -sd_dcells2;
sd_ndcells = -sd_ndcells2;

%% Figure 6N: violin plot of dcells vs non-dcells sensory delay factors
y1 = sd_dcells; 
y2 = sd_ndcells;
hfig = figure('Position',[300 300 300 400]);
hfig.Name = 'Figure 6N: SENSORY DELAY violin plot dist vs non dist';
violin({y1,y2},'edgecolor',lines(1),'facecolor',lines(1),'mc',[],'medc','k');
box off;
xticks([1 2]);
xticklabels({'Dist','Non-Dist'});
set_fig_prefs;
ylabel('Optimal delay factor (sec)');
hold on;
plot(xlim,[0 0],'k:');
legend off;
pval = ranksum(y1,y2);
if pval<0.001
    text(1.5,max(ylim),'***','HorizontalAlignment','center');
elseif pval<0.01
    text(1.5,max(ylim),'**','HorizontalAlignment','center');
elseif pval<0.05
    text(1.5,max(ylim),'*','HorizontalAlignment','center');
else
    text(1.5,max(ylim),'n.s.','HorizontalAlignment','center');
end
plot([1.1 1.9],[0.32 0.32],'k-');

%% Figure 6O: plot by session

mousedate_d = cell(size(sd_uid_dcells));
for cIdx = 1:numel(mousedate_d)
    split_this = strsplit(sd_uid_dcells{cIdx},'_');
    mousedate_d{cIdx} = sprintf('%s_%s',split_this{1},split_this{2});
end
mousedate_nd = cell(size(sd_uid_ndcells));
for cIdx = 1:numel(mousedate_nd)
    split_this = strsplit(sd_uid_ndcells{cIdx},'_');
    mousedate_nd{cIdx} = sprintf('%s_%s',split_this{1},split_this{2});
end

umousedate = intersect(unique(mousedate_d),unique(mousedate_nd));
mouse = cell(size(umousedate));
for sIdx = 1:numel(umousedate)
    split_this = strsplit(umousedate{sIdx},'_');
    mouse{sIdx} = split_this{1};
end

y = nan(numel(umousedate),2);
for sIdx = 1:numel(umousedate)
    keep1 = strcmp(mousedate_d,umousedate{sIdx});
    keep2 = strcmp(mousedate_nd,umousedate{sIdx});
    if sum(keep1)>opt.min_num_dist_cells && sum(keep2)>opt.min_num_dist_cells
        y(sIdx,1) = mean(sd_dcells(keep1));
        y(sIdx,2) = mean(sd_ndcells(keep2));
    end
end
y_diff = y(:,1)-y(:,2);

umousedate = umousedate(~isnan(y_diff));
mouse = mouse(~isnan(y_diff));
y_diff = y_diff(~isnan(y_diff));

% plot
hfig = figure('Position',[300 300 300 400]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6O: SENSORY DELAY dist vs non dist by session';
% myscatter(randn(size(y_diff)),y_diff*1000,lines(1),0.4,60);
umouse_all = unique(dt.Mouse);
plot_col = cbrewer('qual','Set1',numel(umouse_all)); 
for sIdx = 1:numel(y_diff)
    col_idx = find(strcmp(umouse_all,mouse{sIdx}));
    myscatter(randn(1,1),y_diff(sIdx)*1000,plot_col(col_idx,:),0.4,60);
end
xlim([-10 10]);
plot(xlim,[0 0],'k:');
h = gca; h.XAxis.Visible = 'off';
ylabel('Optimal delay (ms), Dist - NonDist');
pval_sesh = signrank(y_diff);
if pval_sesh<0.001
    text(1.5,max(ylim),'***','HorizontalAlignment','center');
elseif pval_sesh<0.01
    text(1.5,max(ylim),'**','HorizontalAlignment','center');
elseif pval_sesh<0.05
    text(1.5,max(ylim),'*','HorizontalAlignment','center');
else
    text(1.5,max(ylim),'n.s.','HorizontalAlignment','center');
end