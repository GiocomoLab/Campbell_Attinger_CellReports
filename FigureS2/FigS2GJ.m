% Makes the following figures:
% Fig S2G,J
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

% K means clustering
rep_cluster = nan(size(score,1),1);
rep_cluster(score(:,1)>opt.score_thresh) = 1;
rep_cluster(score(:,1)<opt.score_thresh) = 2;

%% Figure S2G: plot PC1 score by mouse

br = brain_region_allreps;
mouse = cell(size(rep_id_allreps));
for i = 1:numel(mouse)
    mouse_this = strsplit(rep_id_allreps{i},'_');
    mouse_this = mouse_this{1};
    if contains(mouse_this,'np')
        mouse_this = mouse_this(3:4);
    end
    mouse{i} = mouse_this;
end
mouse_uniq = [unique(mouse(strcmp(br,'MEC'))); unique(mouse(~strcmp(br,'MEC')))];

plot_col = cbrewer('qual','Set2',3,'pchip');
hfig = figure('Position',[500 500 800 400]); hold on;
hfig.Name = 'Figure S2G: PC1 score by mouse';
ctr = 1;
for i = 1:numel(mouse)
    br_idx = find(strcmp(opt.brain_regions,br{i}));
    mouse_idx = find(strcmp(mouse_uniq,mouse{i}));
    col_this = plot_col(br_idx,:);
    scatter(mouse_idx,score(i,1),'MarkerEdgeColor',col_this,'MarkerFaceColor',col_this);
end
xticks(1:numel(mouse_uniq));
xticklabels(mouse_uniq);
xtickangle(90);
xlabel('Mouse');
ylabel('PC1 score');
plot(xlim,[0.4 0.4],'r--');

%% Figure S2J left: avg corrmat and individual session corrmats for MEC baseline data
corrmat_this = rez_all{1}.corrmat_avg_stab_bl;
hfig = figure('Position',[300 300 400 400]); hold on;
hfig.Name = 'Figure S2J left: MEC baseline corrmat';
imagesc(squeeze(nanmean(corrmat_this)));
caxis([0 0.7]);
set(gca,'YDir','normal');
colorbar;

% patches indicating gain value
for tr = 1:num_tr
    patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(1,100),...
        'EdgeColor','none');
    patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(1,100),...
        'EdgeColor','none');
end
xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);

% white lines delineated gain change trials
tr = 0.5+[opt.num_tr_bl];
plot([0.5 num_tr+0.5],[tr tr],'w-');
plot([tr tr],[0.5 num_tr+0.5],'w-');

xlim([-1 num_tr+0.5]);
ylim([-1 num_tr+0.5]);
axis square;
axis off

%% Figure S2J right: individual corrmat plots MEC baseline

fig_counter = 5;
num_rows = 7;

corrmat_this = rez_all{1}.corrmat_avg_stab_bl;
corrmat_this = corrmat_this(:,tf);
corrmat_this = corrmat_this-repmat(mean(corrmat_this),size(corrmat_this,1),1); % center data
score_this = corrmat_this * coeff(:,1); % project onto PC1
[~,sort_idx] = sort(score_this,'descend');
corrmat_sort = rez_all{1}.corrmat_avg_stab_bl(sort_idx,:,:);
    
% figure
num_cols = ceil(size(corrmat_sort,1)/num_rows);
hfig = figure('Position',[100 100 100*num_cols 100*num_rows],'Renderer','Painters');
hfig.Name = 'Figure S2J right: trial corr mat indiv sessions MEC baseline';
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