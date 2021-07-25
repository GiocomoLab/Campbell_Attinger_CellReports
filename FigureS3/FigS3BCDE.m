% Makes the following figures:
% Fig S3B-E
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.sessions = fullfile(paths.intermediate_data,'session_lists');

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% analysis options
opt = load_default_opt;
opt.gain = 0.8;
opt.num_tr_bl = 6;
opt.stab_thresh = 0.5;
opt.min_num_stab_cells = 5;
opt.min_frac_stab_cells = 0; % had 0.2 before, decided to remove this requirement for simplicity
opt.score_thresh = 0.4; % for defining stable vs remapping blocks
opt.brain_regions = {'MEC','VISp','RSP'};
opt.session_files = {'small_gain_MEC','small_gain_VISp','small_gain_RSP'};

%% compute trial corr matrix for each session
% remember that this takes max corr over -opt.max_lag:opt.max_lag
rez_all = compute_trial_corr_mat_all_brain_regions(paths,opt,cell_info);

%% PCA and clustering

% gather data
trgain = [ones(1,opt.num_tr_bl) opt.gain*ones(1,4) ones(1,opt.num_tr_bl)];
num_tr = 2*opt.num_tr_bl + opt.num_tr_gc;
tf = triu(true(num_tr),1);
corrmat_allreps = [];
rep_id_allreps = [];
brain_region_allreps = [];
for i = 1:numel(opt.brain_regions)
    corrmat_allreps = [corrmat_allreps; rez_all{i}.corrmat_avg_stab];
    rep_id_allreps = [rep_id_allreps; rez_all{i}.rep_id_stab];
    num_reps_this = numel(rez_all{i}.rep_id_stab);
    brain_region_allreps = [brain_region_allreps; repmat({rez_all{i}.brain_region},num_reps_this,1)];
end

% PCA
[coeff,score,~,~,expl] = pca(corrmat_allreps(:,tf));

% K means clustering
rep_cluster = nan(size(score,1),1);
rep_cluster(score(:,1)>opt.score_thresh) = 1;
rep_cluster(score(:,1)<opt.score_thresh) = 2;

%% cluster baseline MEC data

% assign MEC baseline reps to clusters
corrmat_this = rez_all{1}.corrmat_avg_stab_bl;
corrmat_this = corrmat_this(:,tf);
corrmat_this = corrmat_this-repmat(mean(corrmat_this),size(corrmat_this,1),1); % center data
pc_proj = corrmat_this * coeff; % project into PCA space
rep_cluster_mecBL = nan(size(pc_proj,1),1);% split into stable and remapping clusters based on threshold
rep_cluster_mecBL(pc_proj(:,1)>opt.score_thresh) = 1; 
rep_cluster_mecBL(pc_proj(:,1)<opt.score_thresh) = 2;
% add to list
brain_region_with_bl = [brain_region_allreps; repmat({'MEC control'},numel(rep_cluster_mecBL),1)];
rep_cluster_with_bl = [rep_cluster; rep_cluster_mecBL];

%% Figure S3B: shiftmat baseline MEC
shiftmat_this = rez_all{1}.shiftmat_avg_stab_bl;
rep_cluster_this = rep_cluster_with_bl(strcmp(brain_region_with_bl,'MEC control'));

hfig_bl = figure; hold on;
hfig_bl.Name = 'Figure S3B: Shiftmat MEC baseline';
avg_this = squeeze(mean(shiftmat_this(rep_cluster_this==1,:,:)));
imagesc(avg_this,'XData',1:num_tr,'YData',1:num_tr);
colorbar
caxis([-5 5]);
xticks([]);
yticks([]);
title('MEC baseline');

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

%% Figure S3C: shift comparisons across reps: PEARSON CORR

hfig = figure('Position',[200 200 1000 400]);
hfig.Name = 'Figure S3C: shift block 1 vs block 2 PEARSON';
set(gcf,'Visible','on');
for br = 1:numel(opt.brain_regions)
    rez_this = rez_all{br};
    rep_id_stab = rep_id_allreps(rep_cluster==1 & ismember(rep_id_allreps,rez_this.uniq_rep_id));
    shiftmat_this = rez_this.shiftmat(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:,:);
    rep_info_this = rez_this.rep_info(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:);
    
    % get unique cell ids
    mouse_session_this = get_mouse_session(rep_info_this.Session); 
    cell_id_this = get_cell_id(rep_info_this.Session,rep_info_this.CellID);  
    shift_this = squeeze(nanmean(shiftmat_this(:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),2));

    % only keep reps with at least two reps from same session
    grp_cell = categorical(cell_id_this);
    grp_rep = categorical(rep_info_this.rep_id); 
    grp_session = categorical(mouse_session_this);
    uniq_sesh = unique(grp_session);
    uniq_rep = [];
    rep_pairs = [];
    for i = 1:numel(uniq_sesh)
        reps = unique(grp_rep(grp_session==uniq_sesh(i)));
        if numel(reps)>1
            uniq_rep = [uniq_rep; reps];
            for j = 1:numel(reps)-1
                for k = j+1:numel(reps)
                    rep_pairs = [rep_pairs; reps(j) reps(k)];
                end
            end
        end
    end    
    keep = ismember(grp_rep,uniq_rep);
    y = shift_this(keep,:);
    y = nanmean(y,2);
    grp_cell = grp_cell(keep);
    grp_rep = grp_rep(keep);
    grp_session = grp_session(keep);    
    
    y2 = [];
    for i = 1:size(rep_pairs,1)
        cells1 = grp_cell(grp_rep==rep_pairs(i,1));
        cells2 = grp_cell(grp_rep==rep_pairs(i,2));
        cells_intersect = intersect(cells1,cells2);
        idx1 = ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,1);
        idx2 = ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,2);
        y2_this = [y(idx1)-nanmean(y(idx1)) y(idx2)-nanmean(y(idx2))];
        y2 = [y2; y2_this];
    end
    
    % plot
    subplot(1,numel(opt.brain_regions),br); hold on;
    scatter(y2(:,1),y2(:,2),'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);    
    r = corrcoef(y2(:,1),y2(:,2));
    lm = fitlm(y2(:,1),y2(:,2));   
    % plot linear fit with confidence intervals
    xpred = linspace(-15,15,100)';
    [ypred, yci] = predict(lm,xpred);
    plot(xpred,ypred,'r-');
    plot(xpred,yci(:,1),'r--');
    plot(xpred,yci(:,2),'r--');
    xlim([-15 15]);
    ylim([-15 15]);
    plot([0 0],ylim(),'k:');
    plot(xlim(),[0 0],'k:');
    xlabel('Map shift, block 1 (mean subtracted)');
    ylabel('Map shift, block 2 (mean subtracted)');
    title(sprintf('%s: n=%d, r=%0.3f, p=%0.1e',opt.brain_regions{br},size(y2,1),r(1,2),lm.coefTest));
    axis square;
end
%% Figure S3D: shift comparisons across reps: NORMALIZED RANK

hfig = figure('Position',[200 200 1000 400]);
hfig.Name = 'Figure S3D: shift block 1 vs block 2 RANK';
set(gcf,'Visible','on');
for br = 1:numel(opt.brain_regions)
    rez_this = rez_all{br};
    rep_id_stab = rep_id_allreps(rep_cluster==1 & ismember(rep_id_allreps,rez_this.uniq_rep_id));
    shiftmat_this = rez_this.shiftmat(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:,:);
    rep_info_this = rez_this.rep_info(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:);
    
    % get unique cell ids
    mouse_session_this = get_mouse_session(rep_info_this.Session); 
    cell_id_this = get_cell_id(rep_info_this.Session,rep_info_this.CellID);
    shift_this = squeeze(nanmean(shiftmat_this(:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),2));
    
    % only keep reps with at least two reps from same session
    grp_cell = categorical(cell_id_this);
    grp_rep = categorical(rep_info_this.rep_id); 
    grp_session = categorical(mouse_session_this);
    uniq_sesh = unique(grp_session);
    uniq_rep = [];
    rep_pairs = [];
    for i = 1:numel(uniq_sesh)
        reps = unique(grp_rep(grp_session==uniq_sesh(i)));
        if numel(reps)>1
            uniq_rep = [uniq_rep; reps];
            for j = 1:numel(reps)-1
                for k = j+1:numel(reps)
                    rep_pairs = [rep_pairs; reps(j) reps(k)];
                end
            end
        end
    end    
    keep = ismember(grp_rep,uniq_rep);
    y = shift_this(keep,:);
    y = nanmean(y,2);
    grp_cell = grp_cell(keep);
    grp_rep = grp_rep(keep);
    grp_session = grp_session(keep); 
    
    % convert shifts to rank order
    y_rank = [];
    for i = 1:size(rep_pairs,1)
        cells1 = grp_cell(grp_rep==rep_pairs(i,1));
        cells2 = grp_cell(grp_rep==rep_pairs(i,2));
        cells_intersect = intersect(cells1,cells2);
        idx1 = find(ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,1));
        idx2 = find(ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,2));
        [~,sort_idx1] = sort(y(idx1));
        [~,sort_idx2] = sort(y(idx2));
        r1 = 1:numel(sort_idx1);
        r2 = 1:numel(sort_idx2);
        r1(sort_idx1) = r1;
        r2(sort_idx2) = r2;
        r1 = r1/numel(r1);
        r2 = r2/numel(r2);
        y_rank = [y_rank; [r1' r2']];
    end        
    
    % plot
    subplot(1,numel(opt.brain_regions),br); hold on;
    scatter(y_rank(:,1),y_rank(:,2),'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);   
    [r,p] = corr(y_rank(:,1),y_rank(:,2),'Type','spearman');
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    xlabel('Normalized rank, block 1');
    ylabel('Normalized rank, block 2');
    title(sprintf('%s: n=%d, r=%0.3f, p=%0.1e',opt.brain_regions{br},size(y_rank,1),r,p));
    axis square;
end

%% Figure S3E: shift comparisons across reps: BAR PLOTS (mean shift change)

for br = 1:numel(opt.brain_regions)
    rez_this = rez_all{br};
    rep_id_stab = rep_id_allreps(rep_cluster==1 & ismember(rep_id_allreps,rez_this.uniq_rep_id));
    shiftmat_this = rez_this.shiftmat(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:,:);
    rep_info_this = rez_this.rep_info(ismember(rez_this.rep_info.rep_id,rep_id_stab) & rez_this.rep_info.Stab_BL>opt.stab_thresh,:);
    
    % get unique cell ids
    mouse_session_this = get_mouse_session(rep_info_this.Session); 
    cell_id_this = get_cell_id(rep_info_this.Session,rep_info_this.CellID);
    shift_this = squeeze(nanmean(shiftmat_this(:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),2));
    
    % only keep reps with at least two reps from same session
    grp_cell = categorical(cell_id_this);
    grp_rep = categorical(rep_info_this.rep_id); 
    grp_session = categorical(mouse_session_this);
    uniq_sesh = unique(grp_session);
    uniq_rep = [];
    rep_pairs = [];
    for i = 1:numel(uniq_sesh)
        reps = unique(grp_rep(grp_session==uniq_sesh(i)));
        if numel(reps)>1
            uniq_rep = [uniq_rep; reps];
            for j = 1:numel(reps)-1
                for k = j+1:numel(reps)
                    rep_pairs = [rep_pairs; reps(j) reps(k)];
                end
            end
        end
    end    
    keep = ismember(grp_rep,uniq_rep);
    y = shift_this(keep,:);
    y = nanmean(y,2);
    grp_cell = grp_cell(keep);
    grp_rep = grp_rep(keep);
    grp_session = grp_session(keep);   
    
    % plot shift diff across blocks as bar plots    
    mean_diff = nan(numel(uniq_sesh),1);
    sem_diff = nan(numel(uniq_sesh),1);
    pval_diff = nan(numel(uniq_sesh),1);
    for i = 1:size(rep_pairs,1)
        cells1 = grp_cell(grp_rep==rep_pairs(i,1));
        cells2 = grp_cell(grp_rep==rep_pairs(i,2));
        cells_intersect = intersect(cells1,cells2);
        idx1 = find(ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,1));
        idx2 = find(ismember(grp_cell,cells_intersect) & grp_rep==rep_pairs(i,2));
        y1 = y(idx1);
        y2 = y(idx2);
        diff_this = y2-y1;
        mean_diff(i) = mean(diff_this);
        sem_diff(i) = std(diff_this)/sqrt(numel(diff_this));
        pval_diff(i) = signrank(y1,y2);
    end        
    
    hfig = figure('Position',[200 200 250 80]); hold on;
    hfig.Name = sprintf('Figure S3E: shifts block1 vs block2 individual sessions BAR %s',opt.brain_regions{br});
    set(gcf,'Visible','on');
    bar(mean_diff,'k');
    errorbar(mean_diff,sem_diff,'k');
    xlim([0 13]);
    xticks([]);
    ylabel(sprintf('Shift (cm)\nblock2-block1'))
    title(opt.brain_regions{br});
    % plot pvals
    axmax = max(ylim);
    for i = 1:numel(pval_diff)
        if pval_diff(i)<0.05
            text(i,axmax+0.3,'*','HorizontalAlignment','center');
        end
    end
        
end