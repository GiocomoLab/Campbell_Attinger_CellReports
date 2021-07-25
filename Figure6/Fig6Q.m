% Makes the following figures:
% Fig 6Q
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
dist_tuning = load(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',binsize,smoothsig)),'dist_tuning','opt');
opt = dist_tuning.opt;
dist_tuning = dist_tuning.dist_tuning;

% analysis options
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence
opt.min_num_cells = 30; % session has to have at least this number of dist-cells AND non-dist-cells to be included here
% downsamples data to have exactly this number of cells in each group for
% each session

% at most one of these can be true:
% (NOTE: Modifying this will give you either Figure 6Q or the controls in Figure S5C)
opt.match_depth = false;
opt.min_depth_diff = 25; % (in microns) minimum distance between cells when matching depth
opt.match_fr = false; 

%% dimensionality of distance vs non-distance tuned neurons: 

%% 1. pick sessions to analyze
dist_cells = dist_tuning.pval<opt.pval_cutoff & dist_tuning.prom>opt.min_prom;
[uniq_sesh,uIdx] = unique(dist_tuning.Session);
mouse = dist_tuning.Mouse(uIdx);
num_cells = nan(numel(uniq_sesh),2); % 1st col: dist-tuned, 2nd col: non-dist-tuned
for i = 1:numel(uniq_sesh)
    num_cells(i,1) = sum(strcmp(dist_tuning.Session,uniq_sesh{i}) & dist_cells);
    num_cells(i,2) = sum(strcmp(dist_tuning.Session,uniq_sesh{i}) & ~dist_cells);
end
keep_sesh = num_cells(:,1)>=opt.min_num_cells & num_cells(:,2)>=opt.min_num_cells; 
sesh_to_analyze = uniq_sesh(keep_sesh);
mouse = mouse(keep_sesh);
num_cells = num_cells(keep_sesh,:);

%% 2. assess dimensionality of dist cells vs non-dist cells using PCA

expl_dist_all = nan(numel(sesh_to_analyze),opt.min_num_cells);
expl_non_dist_all = nan(numel(sesh_to_analyze),opt.min_num_cells);

% variables to control for
mean_fr_dist = nan(numel(sesh_to_analyze),opt.min_num_cells);
mean_fr_non_dist = nan(numel(sesh_to_analyze),opt.min_num_cells);
depth_dist = nan(numel(sesh_to_analyze),opt.min_num_cells);
depth_non_dist = nan(numel(sesh_to_analyze),opt.min_num_cells);
intercell_dist = nan(numel(sesh_to_analyze),opt.min_num_cells,opt.min_num_cells);
intercell_non_dist = nan(numel(sesh_to_analyze),opt.min_num_cells,opt.min_num_cells);

pb = ParforProgressbar(numel(sesh_to_analyze));
parfor i = 1:numel(sesh_to_analyze)
    % load data
    dat = load(fullfile(paths.data,sesh_to_analyze{i}));   
    
    % get dist and non-dist cell ids
    keep_dist = strcmp(dist_tuning.Session,sesh_to_analyze{i}) & dist_cells;
    keep_non_dist = strcmp(dist_tuning.Session,sesh_to_analyze{i}) & ~dist_cells;
    cell_id_dist = dist_tuning.CellID(keep_dist);
    cell_id_non_dist = dist_tuning.CellID(keep_non_dist);   
    depth_dist_this = dist_tuning.DepthAdjusted(keep_dist);
    depth_non_dist_this = dist_tuning.DepthAdjusted(keep_non_dist);
    
    if ~all(isnan(depth_dist_this)) % some sessions have no depth info
    
        % compute fr mat to get mean FR for each cell
        fr_mat_dist = calcFRVsTime(cell_id_dist,dat,opt);
        fr_mat_non_dist = calcFRVsTime(cell_id_non_dist,dat,opt); 
        mean_fr_dist_this = mean(fr_mat_dist,2);
        mean_fr_non_dist_this = mean(fr_mat_non_dist,2);

        % find closest match in depth OR firing ratefor each neuron  
        cell_id = {cell_id_dist, cell_id_non_dist};
        if opt.match_depth           
            value_to_match = {depth_dist_this, depth_non_dist_this};  
        elseif opt.match_fr
            value_to_match = {mean_fr_dist_this, mean_fr_non_dist_this};  
        end
        if opt.match_depth || opt.match_fr
            keep = cell(1,2);
            if numel(cell_id_dist)<=numel(cell_id_non_dist)
                grp1 = 1;
                grp2 = 2;
            else
                grp1 = 2;
                grp2 = 1;
            end
            keep{grp2} = nan(size(cell_id{grp1}));
            diff_all = nan(size(cell_id{grp1}));
            for j = 1:numel(cell_id{grp1})
                diff_this = abs(value_to_match{grp2}-value_to_match{grp1}(j));
                if opt.match_depth
                    diff_this(diff_this<=opt.min_depth_diff)=Inf;
                end
                [min_diff,min_idx] = min(diff_this);
                keep{grp2}(j) = min_idx;
                diff_all(j) = min_diff;
                value_to_match{grp2}(min_idx) = nan; % remove this cell from the pool
            end
            [~,sort_idx] = sort(diff_all);           
            % take the top opt.min_num_cells closests matching cell pairs
            keep{grp1} = sort_idx(1:opt.min_num_cells);
            keep{grp2} = keep{grp2}(keep{grp1});   
            keep_dist = keep{1};
            keep_non_dist = keep{2};
        else
            % no matching
            % randomly downsample to have opt.min_num_cells in each group
            [~,keep_dist] = datasample(cell_id_dist,opt.min_num_cells,'Replace',false);
            [~,keep_non_dist] = datasample(cell_id_non_dist,opt.min_num_cells,'Replace',false);      
        end

        % keep only selected cells
        cell_id_dist= cell_id_dist(keep_dist);
        cell_id_non_dist = cell_id_non_dist(keep_non_dist);          
        fr_mat_dist = fr_mat_dist(keep_dist,:);
        fr_mat_non_dist = fr_mat_non_dist(keep_non_dist,:);
        depth_dist_this = depth_dist_this(keep_dist);
        depth_non_dist_this = depth_non_dist_this(keep_non_dist);  
        mean_fr_dist_this = mean(fr_mat_dist,2);
        mean_fr_non_dist_this = mean(fr_mat_non_dist,2);

        if numel(cell_id_dist)>=opt.min_num_cells && numel(cell_id_non_dist)>=opt.min_num_cells

            % take zscore
            fr_mat_dist_zscore = my_zscore(fr_mat_dist);
            fr_mat_non_dist_zscore = my_zscore(fr_mat_non_dist);

            % pca on zscored firing rate matrix
            [~,~,~,~,expl_dist] = pca(fr_mat_dist_zscore');
            [~,~,~,~,expl_non_dist] = pca(fr_mat_non_dist_zscore');

            % save variance explained
            expl_dist_all(i,:) = expl_dist;
            expl_non_dist_all(i,:) = expl_non_dist;

            % save nuisance variables (to check later)
            mean_fr_dist(i,:) = mean_fr_dist_this;
            mean_fr_non_dist(i,:) = mean_fr_non_dist_this;
            depth_dist(i,:) = depth_dist_this;
            depth_non_dist(i,:) = depth_non_dist_this;
            intercell_dist_this = nan(opt.min_num_cells,opt.min_num_cells);
            intercell_non_dist_this = nan(opt.min_num_cells,opt.min_num_cells);
            for j = 1:opt.min_num_cells-1
                for k = j+1:opt.min_num_cells
                    intercell_dist_this(j,k) = abs(depth_dist_this(j)-depth_dist_this(k));
                    intercell_non_dist_this(j,k) = abs(depth_non_dist_this(j)-depth_non_dist_this(k));
                end
            end
            intercell_dist(i,:,:) = intercell_dist_this;
            intercell_non_dist(i,:,:) = intercell_non_dist_this;

        end
    end

    pb.increment();
end

%% reshape variables

mean_fr = nan(numel(sesh_to_analyze)*opt.min_num_cells,2);
depth = nan(numel(sesh_to_analyze)*opt.min_num_cells,2);
intercell = nan(numel(sesh_to_analyze)*opt.min_num_cells*opt.min_num_cells,2);

mean_fr(:,1) = reshape(mean_fr_dist,numel(mean_fr_dist),1);
mean_fr(:,2) = reshape(mean_fr_non_dist,numel(mean_fr_non_dist),1);
mean_fr = mean_fr(all(~isnan(mean_fr),2),:);
depth(:,1) = reshape(depth_dist,numel(depth_dist),1);
depth(:,2) = reshape(depth_non_dist,numel(depth_non_dist),1);
depth = depth(all(~isnan(depth),2),:);
intercell(:,1) = reshape(intercell_dist,numel(intercell_dist),1);
intercell(:,2) = reshape(intercell_non_dist,numel(intercell_non_dist),1);
intercell = intercell(all(~isnan(intercell),2),:);

%% 3. plot results
hfig = figure('Position',[400 400 300 300]); hold on;
if opt.match_depth
    hfig.Name = 'Figure S5C middle: Dimensionality dist vs non-dist_match depth';
elseif opt.match_fr
    hfig.Name = 'Figure S5C bottom: Dimensionality dist vs non-dist_match fr';
else
    hfig.Name = 'Figure 6Q and Figure S5C top: Dimensionality dist vs non-dist_no match';
end


% cumulative sum of % expl
y1 = cumsum(expl_dist_all,2);
y2 = cumsum(expl_non_dist_all,2);
keep = all(~isnan(y2),2) & all(~isnan(y1),2);
y1 = y1(keep,:);
y2 = y2(keep,:);
sesh_keep = sesh_to_analyze(keep);
mouse_keep = mouse(keep);

% distance cells
mu1 = mean(y1);
sem1 = std(y1)/sqrt(size(y1,1));
errorbar(mu1,sem1,'b');

% non-distance cells
mu2 = mean(y2);
sem2 = std(y2)/sqrt(size(y2,1));
errorbar(mu2,sem2,'k');

% plot significance
diff = y2-y1;
for i = 1:size(diff,2)
    p = signrank(diff(:,i));
    if p<0.001
        plot([i i i],[98 99 100],'k.');
    elseif p<0.01
        plot([i i],[99 100],'k.');
    elseif p<0.05
        plot(i,100,'k.');
    end
end

legend({'Dist-tuned','Non-dist-tuned'},'location','southeast');
xlabel('Num. PCs');
ylabel('% Var. Explained');
title(sprintf('n = %d sessions from %d mice',numel(sesh_keep),numel(unique(mouse_keep))));
set_fig_prefs;