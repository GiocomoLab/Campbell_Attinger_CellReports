% Makes the following figures:
% Fig 2D-G
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

%restoredefaultpath;

% add helper functions to path
%addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.sessions = fullfile(paths.intermediate_data,'session_lists');

% load cell_info table
load(fullfile(paths.intermediate_data,fullfile('cell_info','cell_info_Apr2020'))); 

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

% Classify reps as stable (cluster 1) or remap (cluster 2)
rep_cluster = nan(size(score,1),1);
rep_cluster(score(:,1)>opt.score_thresh) = 1;
rep_cluster(score(:,1)<opt.score_thresh) = 2;

%% Figure 2D: Avg corrmat for each brain region

hfig = figure('Position',[300 300 1000 400]);
hfig.Name = 'Figure 2D: avg corrmats by brain region (MEC, RSC, V1)';
for br = 1:numel(opt.brain_regions)
    corrmat_this = rez_all{br}.corrmat_avg_stab;
    
    subplot(1,numel(opt.brain_regions),br); hold on;
    avg_this = squeeze(mean(corrmat_this));
    imagesc(avg_this);
    colorbar
    caxis([0 0.7]);
    xticks([]);
    yticks([]);
    
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
end

%% Figure 2E: Indiv corrmat plot
trgain = [ones(1,opt.num_tr_bl) opt.gain*ones(1,4) ones(1,opt.num_tr_bl)];
fig_counter = 5;
num_rows = 7;
for br = 1:numel(opt.brain_regions)
    corrmat_sort = rez_all{br}.corrmat_avg_stab;
    score_this = score(strcmp(brain_region_allreps,opt.brain_regions{br}),1);
    [~,sort_idx] = sort(score_this,'descend');
    corrmat_sort = corrmat_sort(sort_idx,:,:);
    
    % figure
    num_cols = ceil(size(corrmat_sort,1)/num_rows);
    hfig = figure('Position',[100 100 100*num_cols 100*num_rows],'Renderer','Painters');
    hfig.Name = sprintf('Figure 2E: trial corr mat indiv sessions %s',rez_all{br}.brain_region);
    fig_counter = fig_counter+1;
    ha = tight_subplot(num_rows,num_cols,[0.01 0.01],[0.01 0.01],[0.01 0.01]);
    for i = 1:size(corrmat_sort,1)
        axes(ha(i)); hold on;
        imagesc(squeeze(corrmat_sort(i,:,:)));
        caxis([0 0.7]);
        xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);
        axis square;
        axis off;
    end
end

%% Figure 2F: histograms of PC1 score

hfig_histograms = figure('Position',[300 300 625 350]);
hfig_histograms.Name = 'Figure 2F: histograms of PC1 scores by brian region';
plot_col = cbrewer('qual','Set2',3,'pchip');
subplot(1,numel(opt.brain_regions)+2,1); hold on;
histogram(score(:,1),-3:0.2:4,'FaceColor','k','Orientation','horizontal');
plot(xlim,[0.4 0.4],'r-');
grid on;
box off;
ylabel('PC1 score');
title('All');
set(gca,'FontSize',12);
for i = 1:numel(opt.brain_regions)
    subplot(1,numel(opt.brain_regions)+2,i+1); hold on;
    histogram(score(strcmp(brain_region_allreps,opt.brain_regions{i}),1),-3:0.2:4,'FaceColor',plot_col(i,:),'Orientation','horizontal');
    plot(xlim,[opt.score_thresh opt.score_thresh],'r-');
    yticklabels('');
    grid on;
    box off;
    title(opt.brain_regions{i});
    set(gca,'FontSize',12);
end

% MEC baseline score
corrmat_this = rez_all{1}.corrmat_avg_stab_bl;
corrmat_this = corrmat_this(:,tf);
corrmat_this = corrmat_this-repmat(mean(corrmat_this),size(corrmat_this,1),1); % center data
pc_proj = corrmat_this * coeff(:,1); % project onto PC1

subplot(1,numel(opt.brain_regions)+2,numel(opt.brain_regions)+2); hold on;
histogram(pc_proj,-3:0.2:4,'FaceColor','k','Orientation','horizontal');
plot(xlim,[0.4 0.4],'r-');
yticklabels('');
grid on;
box off;
title('MECb');
set(gca,'FontSize',12);

%% Figure 2G top: Avg corrmat of each cluster

hfig = figure('Position',[300 300 400 150]);
hfig.Name = 'Figure 2G top: avg cluster corrmats';
for i = 1:2
    subplot(1,2,i); hold on;
    avg_this = squeeze(mean(corrmat_allreps(rep_cluster==i,:,:)));
    imagesc(avg_this);
    caxis([0 0.7]);
    xticks([]);
    yticks([]);
    
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
end

%% Figure 2G bottom: Bar plots of number of reps in each cluster

% assign MEC baseline reps to clusters
corrmat_this = rez_all{1}.corrmat_avg_stab_bl;
corrmat_this = corrmat_this(:,tf);
corrmat_this = corrmat_this-repmat(mean(corrmat_this),size(corrmat_this,1),1); % center data
pc_proj = corrmat_this * coeff; % project into PCA space
rep_cluster_mecBL = nan(size(pc_proj,1),1);% split into stable and remapping clusters based on threshold
rep_cluster_mecBL(pc_proj(:,1)>opt.score_thresh) = 1; 
rep_cluster_mecBL(pc_proj(:,1)<opt.score_thresh) = 2;
% add to list
brain_region_this = [brain_region_allreps; repmat({'MEC control'},numel(rep_cluster_mecBL),1)];
rep_cluster_this = [rep_cluster; rep_cluster_mecBL];

% plot
hfig = figure('Position',[400 400 500 250]);
hfig.Name = 'Figure 2G bottom: bar graph of num sessions of each type';
bar_data  = crosstab(brain_region_this,rep_cluster_this);
% bar_data_norm = 100*bar_data./(repmat(sum(bar_data,2),1,opt.num_clust));
barh(flipud(bar_data),'stacked');
xlabel('Num. reps');
yticklabels(fliplr({'MEC','V1','RSC','MEC baseline'}));