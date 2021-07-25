% Makes the following figures:
% Fig 3A,C,D
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

%restoredefaultpathe ;

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

%% Figure 3A: spatial cells per block bar plot

num_cells = [];
num_cells_br = []; % brain region
for i = 1:numel(opt.brain_regions)
    rep_info_this = rez_all{i}.rep_info; 
    uniq_rep_this = rez_all{i}.uniq_rep_id;
    num_reps_this = numel(uniq_rep_this);
    num_cells_this = nan(num_reps_this,2);
    for j = 1:num_reps_this        
        num_cells_this(j,1) = sum(strcmp(rep_info_this.rep_id,uniq_rep_this{j}) & ...
            rep_info_this.Stab_BL > opt.stab_thresh);
        num_cells_this(j,2) = sum(strcmp(rep_info_this.rep_id,uniq_rep_this{j}) & ...
            (rep_info_this.Stab_BL <= opt.stab_thresh | isnan(rep_info_this.Stab_BL)));
    end
    [~,sort_idx] = sort(num_cells_this(:,1),'descend');
    num_cells = [num_cells; num_cells_this(sort_idx,:)];
    num_cells_br = [num_cells_br; i*ones(size(num_cells_this,1),1)];
end
num_cells_filt = num_cells(num_cells(:,1)>=opt.min_num_stab_cells,:);
num_cells_br_filt = num_cells_br(num_cells(:,1)>=opt.min_num_stab_cells);

% only keep stable blocks
num_cells_stab = num_cells_filt(rep_cluster==1,:);
num_cells_br_stab = num_cells_br_filt(rep_cluster==1);

num_spatial_cells_plot = cell(numel(opt.brain_regions),1);
means = nan(numel(opt.brain_regions),1);
for i = 1:numel(opt.brain_regions)
    num_cells_this = num_cells_stab(num_cells_br_stab==i,1);
    num_spatial_cells_plot{i} = num_cells_this;
    means(i) = mean(num_cells_this);
    errs(i) = std(num_cells_this)/sqrt(numel(num_cells_this));
end

% plot
hfig = figure('Position',[300 300 300 300]);
hfig.Name = 'Figure 3A: num spatial cells bar plot';
set(gcf,'Visible','on');
bar(means); hold on;
errorbar(means,errs,'k.');
plotSpread(num_spatial_cells_plot);
xticklabels(opt.brain_regions);
ylabel('Num. spatial cells');

%% Figure 3C: Shiftmat plot

hfig = figure('Position',[300 300 1000 400]);
hfig.Name = 'Figure 3C: avg shiftmats by brain region';
set(gcf,'Visible','on');
for br = 1:numel(opt.brain_regions)
    shiftmat_this = rez_all{br}.shiftmat_avg_stab;
    rep_cluster_this = rep_cluster(strcmp(brain_region_allreps,opt.brain_regions{br}));
    
    subplot(1,numel(opt.brain_regions),br); hold on;
    avg_this = squeeze(mean(shiftmat_this(rep_cluster_this==1,:,:)));
    imagesc(avg_this,'XData',1:num_tr,'YData',1:num_tr);
    colorbar
    caxis([-5 5]);
    xticks([]);
    yticks([]);
    title(sprintf('%d stable blocks',sum(rep_cluster_this==1)));

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

%% Figure 3D: Shift comparison bar plot

% gather data
shift_all = [];
shift_brain_region = [];
shift_plot = {};
means = nan(numel(opt.brain_regions),1);
errs = nan(numel(opt.brain_regions),1);
for br = 1:numel(opt.brain_regions)
    shiftmat_this = rez_all{br}.shiftmat_avg_stab;
    rep_cluster_this = rep_cluster(strcmp(brain_region_allreps,opt.brain_regions{br}));
    avg_shift_this = mean(mean(shiftmat_this(rep_cluster_this==1,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),3),2);    
    shift_all = [shift_all; avg_shift_this];
    shift_brain_region = [shift_brain_region; br*ones(size(avg_shift_this,1),1)];
    shift_plot{br} = avg_shift_this;
    means(br) = mean(avg_shift_this);
    errs(br) = std(avg_shift_this)/sqrt(numel(avg_shift_this));
end

% plot
hfig = figure('Position',[300 300 300 300]);
hfig.Name = 'Figure 3D: shift brain region comparison';
set(gcf,'Visible','on');
bar(means); hold on;
errorbar(means,errs,'k.');
plotSpread(shift_plot);
xticklabels(opt.brain_regions);
ylabel('Map shift (cm)');

% anova
anovap = anovan(shift_all,shift_brain_region,'display','off');
if anovap<0.001
    title('***');
elseif anovap<0.01
    title('**');
elseif anovap<0.05
    title('*');
end

