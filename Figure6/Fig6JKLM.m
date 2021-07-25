% Makes the following figures:
% Fig 6J,K,L,M
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = fullfile(paths.intermediate_data,'dark_distance');
paths.sessions = fullfile(paths.intermediate_data,'session_lists');

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% load distance tuning results
dtab = load(fullfile(paths.rez,'dist_tuning_autocorr_bin2_smooth4.mat'));
opt = dtab.opt;
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1;
dtab = dtab.dist_tuning;
dsig = dtab.pval<opt.pval_cutoff & dtab.prom>opt.min_prom;
dcells = dtab.UniqueID(dsig);
opt.min_num_dist_tuned_cells = 5; % for analyses by rep

% load rep clustering results
rc = load(fullfile(paths.intermediate_data,'rep_clusters','rep_clusters.mat')); 
rep_stable = rc.rep_id(rc.rep_cluster==1 & strcmp(rc.brain_region,'MEC'));

% for shift calculations
session_file = 'small_gain_MEC';
opt.dark = false;
opt.gain = 0.8;
opt.max_lag = 30;
opt.num_tr_bl = 6;
opt.stab_thresh = 0.5;
opt.min_num_stab_cells = 5;
opt.min_frac_stab_cells = 0; % had 0.2 before, decided to remove this requirement for simplicity

%% compute trial corr matrix for each session
% remember that this takes max corr over -opt.max_lag:opt.max_lag

rez = struct;
rez.brain_region = opt.brain_region;
rez.session_name = read_session_file(session_file,paths);
[rez.corrmat,rez.shiftmat,rez.rep_info,...
    rez.corrmat_bl,rez.shiftmat_bl,rez.rep_info_bl] = ...
    compute_trial_corr_mat_indiv_neurons_with_baseline(rez.session_name,cell_info,paths,opt);

% GAIN CHANGE
% average by rep 
rez.uniq_rep_id = unique(rez.rep_info.rep_id);
rez.corrmat_avg = nan(numel(rez.uniq_rep_id),size(rez.corrmat,2),size(rez.corrmat,2));
rez.shiftmat_avg = nan(numel(rez.uniq_rep_id),size(rez.corrmat,2),size(rez.corrmat,2));
stab_reps = false(size(rez.corrmat_avg,1),1);
for i = 1:numel(rez.uniq_rep_id)
    keep_cell1 = strcmp(rez.rep_info.rep_id,rez.uniq_rep_id{i});
    keep_cell2 = rez.rep_info.Stab_BL>opt.stab_thresh;
    keep_cell = keep_cell1 & keep_cell2;
    rez.corrmat_avg(i,:,:) = squeeze(nanmean(rez.corrmat(keep_cell,:,:),1));
    rez.shiftmat_avg(i,:,:) = squeeze(nanmean(rez.shiftmat(keep_cell,:,:),1));
    stab_reps(i) = sum(keep_cell)>=opt.min_num_stab_cells && sum(keep_cell)/sum(keep_cell1)>=opt.min_frac_stab_cells;
end
% only keep stable reps
rez.rep_id_stab = rez.uniq_rep_id(stab_reps);
rez.corrmat_avg_stab = rez.corrmat_avg(stab_reps,:,:);
rez.shiftmat_avg_stab = rez.shiftmat_avg(stab_reps,:,:);

% BASELINE
% average by rep 
rez.uniq_rep_id_bl = unique(rez.rep_info_bl.rep_id);
rez.corrmat_avg_bl = nan(numel(rez.uniq_rep_id_bl),size(rez.corrmat_bl,2),size(rez.corrmat_bl,2));
rez.shiftmat_avg_bl = nan(numel(rez.uniq_rep_id_bl),size(rez.corrmat_bl,2),size(rez.corrmat_bl,2));
stab_reps = false(size(rez.corrmat_avg_bl,1),1);
for i = 1:numel(rez.uniq_rep_id_bl)
    keep_cell1 = strcmp(rez.rep_info_bl.rep_id,rez.uniq_rep_id_bl{i});
    keep_cell2 = rez.rep_info_bl.Stab_BL>opt.stab_thresh;
    keep_cell = keep_cell1 & keep_cell2;
    rez.corrmat_avg_bl(i,:,:) = squeeze(nanmean(rez.corrmat_bl(keep_cell,:,:),1));
    rez.shiftmat_avg_bl(i,:,:) = squeeze(nanmean(rez.shiftmat_bl(keep_cell,:,:),1));
    stab_reps(i) = sum(keep_cell)>=opt.min_num_stab_cells && sum(keep_cell)/sum(keep_cell1)>=opt.min_frac_stab_cells;
end
% only keep stable reps
rez.rep_id_stab_bl = rez.uniq_rep_id_bl(stab_reps);
rez.corrmat_avg_stab_bl = rez.corrmat_avg_bl(stab_reps,:,:);
rez.shiftmat_avg_stab_bl = rez.shiftmat_avg_bl(stab_reps,:,:);

%% Figure 6J: plot by cell

% extract average shifts for cells with dark distance data
keep = ismember(rez.rep_info.rep_id,rep_stable) & rez.rep_info.Stab_BL>opt.stab_thresh;
corrmat = rez.corrmat(keep,:,:);
shiftmat = rez.shiftmat(keep,:,:);
rep_info = rez.rep_info(keep,:);
uid = cell(size(rep_info,1),1);
for cIdx = 1:numel(uid)
    split_this = strsplit(rep_info.Session{cIdx},'_');
    uid{cIdx} = sprintf('%s_%s_c%d',split_this{1},split_this{2},rep_info.CellID(cIdx));
end
rep_info.UniqueID = uid;

% only keep those with dark distance data
keep = ismember(rep_info.UniqueID,dtab.UniqueID);
corrmat = corrmat(keep,:,:);
shiftmat = shiftmat(keep,:,:);
rep_info = rep_info(keep,:);
dist_tuned = ismember(rep_info.UniqueID,dcells);

% avg shift
avg_shift = nanmean(nanmean(shiftmat(:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),3),2);

hfig = figure('Position',[300 300 300 400]); 
hfig.Name = 'Figure 6J: shift dist vs non dist Cluster1 by cell';

dcells_this = unique(rep_info.UniqueID(dist_tuned));
ndcells_this = unique(rep_info.UniqueID(~dist_tuned));

y1 = nan(numel(dcells_this),1);
y2 = nan(numel(ndcells_this),1);
for i = 1:numel(dcells_this)
    y1(i) = nanmean(avg_shift(strcmp(rep_info.UniqueID,dcells_this{i})));
end
for i = 1:numel(ndcells_this)
    y2(i) = nanmean(avg_shift(strcmp(rep_info.UniqueID,ndcells_this{i})));
end
y1 = y1(~isnan(y1));
y2 = y2(~isnan(y2));

violin({y1,y2},'mc',[]);
legend off
ylabel('Shift (cm)');
plot(xlim,[0 0],'k--');
xticks([1 2]);
yticks(-30:10:30);
xticklabels({'Dist','Non-Dist'});
box off

pval = ranksum(y1,y2);
plot([1.1 1.9],[max(ylim) max(ylim)],'k');
if pval<0.001
    text(1.5,max(ylim)+2,'***','HorizontalAlignment','center');
elseif pval<0.01
    text(1.5,max(ylim)+2,'**','HorizontalAlignment','center');
elseif pval<0.05
    text(1.5,max(ylim)+2,'*','HorizontalAlignment','center');
else
    text(1.5,max(ylim)+2,'n.s.','HorizontalAlignment','center');
end
set_fig_prefs;

%% Figure 6K: plot by rep
rep_id = unique(rep_info.rep_id);
ndist = nan(numel(rep_id),1);
for rIdx = 1:numel(rep_id)
    ndist(rIdx) = sum(ismember(rep_info.UniqueID,dcells) & strcmp(rep_info.rep_id,rep_id{rIdx}));
end
rep_id = rep_id(ndist>=opt.min_num_dist_tuned_cells);
ndist = ndist(ndist>=opt.min_num_dist_tuned_cells);

mouse = cell(size(rep_id));
for rIdx = 1:numel(rep_id)
    split_this = strsplit(rep_id{rIdx},'_');
    mouse{rIdx} = split_this{1};
end

shift_diff_by_rep = nan(numel(rep_id),1);
pval_by_rep = nan(numel(rep_id),1);
for rIdx = 1:numel(rep_id)
    y1 = avg_shift(dist_tuned & strcmp(rep_info.rep_id,rep_id{rIdx}));
    y2 = avg_shift(~dist_tuned & strcmp(rep_info.rep_id,rep_id{rIdx}));
    shift_diff_by_rep(rIdx) = nanmean(y1)-nanmean(y2);
    pval_by_rep(rIdx) = ranksum(y1,y2);
end

hfig = figure('Position',[400 400 150 300]); hold on; set_fig_prefs;
hfig.Name = 'Figure 6K: shift dist vs non dist Cluster1 by rep';


umouse_all = unique(dtab.Mouse);
plot_col = cbrewer('qual','Set1',numel(umouse_all)); 
for sIdx = 1:numel(shift_diff_by_rep)
    col_idx = find(strcmp(umouse_all,mouse{sIdx}));
    myscatter(randn(1,1),shift_diff_by_rep(sIdx),plot_col(col_idx,:),0.4,60);
end

% keep = pval_by_rep>0.05; N = sum(keep);
% myscatter(randn(N,1),shift_diff_by_rep(keep),'k',0.1,50);
% keep = pval_by_rep<=0.05; N = sum(keep);
% myscatter(randn(N,1),shift_diff_by_rep(keep),'k',0.5,50);
xlim([-20 20]);
plot(xlim,[0 0],'k:');
h = gca; h.XAxis.Visible = 'off';
ylabel('Map shift, Dist - NonDist (cm)');
%legend({'p>0.05','p<0.05'});

%% Figure 6L: plot average corrmats dist vs non-dist

keep = ~ismember(rez.rep_info.rep_id,rep_stable) & rez.rep_info.Stab_BL>opt.stab_thresh;
corrmat = rez.corrmat(keep,:,:);
shiftmat = rez.shiftmat(keep,:,:);
rep_info = rez.rep_info(keep,:);
uid = cell(size(rep_info,1),1);
for cIdx = 1:numel(uid)
    split_this = strsplit(rep_info.Session{cIdx},'_');
    uid{cIdx} = sprintf('%s_%s_c%d',split_this{1},split_this{2},rep_info.CellID(cIdx));
end
rep_info.UniqueID = uid;

% only keep those with dark distance data
keep = ismember(rep_info.UniqueID,dtab.UniqueID);
corrmat = corrmat(keep,:,:);
shiftmat = shiftmat(keep,:,:);
rep_info = rep_info(keep,:);
dist_tuned = ismember(rep_info.UniqueID,dcells);

dcells_this = unique(rep_info.UniqueID(dist_tuned));
ndcells_this = unique(rep_info.UniqueID(~dist_tuned));

mat = cell(2,1);
mat{1} = nan(numel(dcells_this),size(corrmat,2),size(corrmat,3));
mat{2} = nan(numel(ndcells_this),size(corrmat,2),size(corrmat,3));
for i = 1:numel(dcells_this)
    mat{1}(i,:,:) = nanmean(corrmat(strcmp(rep_info.UniqueID,dcells_this{i}),:,:),1);
end
for i = 1:numel(ndcells_this)
    mat{2}(i,:,:) = nanmean(corrmat(strcmp(rep_info.UniqueID,ndcells_this{i}),:,:),1);
end

hfig = figure('Position',[500 500 600 300]);
hfig.Name = 'Figure 6L: corrmats dist vs non dist by cell';
num_tr = 2*opt.num_tr_bl + opt.num_tr_gc;
trgain = [ones(1,opt.num_tr_bl) opt.gain*ones(1,opt.num_tr_gc) ones(1,opt.num_tr_bl)];
titles = {'Dist','Non-Dist'};
for pIdx = 1:2
    subplot(1,2,pIdx); hold on;
    imagesc(squeeze(nanmean(mat{pIdx})));
    caxis([0.2 0.8]);
    axis square
    colorbar
    
    % patches indicating gain value
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor','none');
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor','none');
    end
    xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);
    
    % white lines delineated gain change trials
    for tr = 0.5+[opt.num_tr_bl opt.num_tr_bl+opt.num_tr_gc]
        plot([0.5 num_tr+0.5],[tr tr],'w-');
        plot([tr tr],[0.5 num_tr+0.5],'w-');
    end
    
    xlim([-1 num_tr+0.5]);
    ylim([-1 num_tr+0.5]);
    axis square;
    axis off
    title(sprintf('%s (%d units)',titles{pIdx},size(mat{pIdx},1)));
end

%% Figure 6M: plot avg corr to baseline with error bars

corr_to_bl_dist = nanmean(mat{1}(:,:,1:opt.num_tr_bl),3);
corr_to_bl_ndist = nanmean(mat{2}(:,:,1:opt.num_tr_bl),3);

means = nan(2,size(corr_to_bl_dist,2));
sems = nan(2,size(corr_to_bl_dist,2));
pvals = nan(1,size(corr_to_bl_dist,2));
for trIdx = 1:size(corr_to_bl_dist,2)
    y1 = corr_to_bl_dist(:,trIdx);
    y2 = corr_to_bl_ndist(:,trIdx);
    means(1,trIdx) = nanmean(y1);
    means(2,trIdx) = nanmean(y2);
    sems(1,trIdx) = nanstd(y1)/sqrt(sum(~isnan(y1)));
    sems(2,trIdx) = nanstd(y2)/sqrt(sum(~isnan(y2)));
    pvals(trIdx) = ranksum(y1,y2);
end

hfig = figure('Position',[300 300 400 350]); hold on;
hfig.Name = 'Figure 6M: corr to bl dist vs non dist line with error bar';
plot_col = cbrewer('qual','Set1',3);
for idx = 1:2
    errorbar(1:size(corr_to_bl_dist,2),means(idx,:),sems(idx,:),'Color',plot_col(idx,:));
end
legend({'Dist','Non-Dist'});
xticks([1 16]);
plot([6.5 6.5],ylim,'k--','HandleVisibility','off');
plot([10.5 10.5],ylim,'k--','HandleVisibility','off');
ylim([0.2 0.8]);
xlabel('Trial');
ylabel('Corr. to baseline');

% significance dots
for i = 1:numel(pvals)
    if pvals(i)<0.001
        plot(i,max(ylim),'k.','HandleVisibility','off');
        plot(i,max(ylim)-0.01,'k.','HandleVisibility','off');
        plot(i,max(ylim)-0.02,'k.','HandleVisibility','off');
    elseif pvals(i)<0.01
        plot(i,max(ylim),'k.','HandleVisibility','off');
        plot(i,max(ylim)-0.01,'k.','HandleVisibility','off');
    elseif pvals(i)<0.05
        plot(i,max(ylim),'k.','HandleVisibility','off');
    end
end