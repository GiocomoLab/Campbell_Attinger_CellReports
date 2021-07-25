% Makes the following figures:
% Fig 5F
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% MGC 7/19/2021

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
opt.contr_hi = 100;
opt.contr_lo = 10;
opt.num_tr_bl = 6;
opt.stab_thresh = 0.5;
opt.remap_thresh = 0.5;
opt.max_lag_peak_fr = 30;

%% Compute correlations

% MEC
opt.brain_region = 'MEC';
rez_mec = struct; % struct for holding MEC results
rez_mec.brain_region = opt.brain_region;
rez_mec.session_file = 'small_gain_contrast10_MEC';
rez_mec.session_name = read_session_file(rez_mec.session_file,paths);
[rez_mec.corrmat_hi, rez_mec.corrmat_lo, rez_mec.shiftmat_hi, rez_mec.shiftmat_lo, rez_mec.repinfo_hi, rez_mec.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(rez_mec.session_file,cell_info,opt,paths);

% VISp
opt.brain_region = 'VISp';
rez_vis = struct; % struct for holding MEC results
rez_vis.brain_region = opt.brain_region;
rez_vis.session_file = 'small_gain_contrast10_VISp';
rez_vis.session_name = read_session_file(rez_vis.session_file,paths);
[rez_vis.corrmat_hi, rez_vis.corrmat_lo, rez_vis.shiftmat_hi, rez_vis.shiftmat_lo, rez_vis.repinfo_hi, rez_vis.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(rez_vis.session_file,cell_info,opt,paths);

% RSP
opt.brain_region = 'RSP';
rez_rsp = struct; % struct for holding MEC results
rez_rsp.brain_region = opt.brain_region;
rez_rsp.session_file = 'small_gain_contrast10_RSP';
rez_rsp.session_name = read_session_file(rez_rsp.session_file,paths);
[rez_rsp.corrmat_hi, rez_rsp.corrmat_lo, rez_rsp.shiftmat_hi, rez_rsp.shiftmat_lo, rez_rsp.repinfo_hi, rez_rsp.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(rez_rsp.session_file,cell_info,opt,paths);


%% find cells that pass criteria

% MEC

% dim1=sesh num, 2=cell id, 3=rep num, 4=gc trial num
rez_mec.num_cells_total = 0;
rez_mec.cell_id_hi = [];
rez_mec.cell_id_lo = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_mec.repinfo_hi)    
    if ~isempty(rez_mec.corrmat_hi{i})

        cell_id_this = rez_mec.repinfo_hi{i}.CellID{1};
        rez_mec.num_cells_total = rez_mec.num_cells_total + numel(cell_id_this);
        cell_id_hi_this = [];
        cell_id_lo_this = [];
        
        % high contrast        
        trial_corr_mat_this = rez_mec.corrmat_hi{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_mec.repinfo_hi{i},1)
                stab = rez_mec.repinfo_hi{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_hi_this = [cell_id_hi_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end       

        % low contrast
        trial_corr_mat_this = rez_mec.corrmat_lo{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_mec.repinfo_lo{i},1)
                stab = rez_mec.repinfo_lo{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_lo_this = [cell_id_lo_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end

        % add to list
        rez_mec.cell_id_hi = [rez_mec.cell_id_hi; cell_id_hi_this];
        rez_mec.cell_id_lo = [rez_mec.cell_id_lo; cell_id_lo_this];
    end
end
% create id strings for cells
rez_mec.cell_name_hi = [];
rez_mec.cell_name_lo = [];
for i = 1:size(rez_mec.cell_id_hi,1)
    rez_mec.cell_name_hi = [rez_mec.cell_name_hi; {sprintf('%d_%d',rez_mec.cell_id_hi(i,1),rez_mec.cell_id_hi(i,2))}];
end
for i = 1:size(rez_mec.cell_id_lo,1)
    rez_mec.cell_name_lo = [rez_mec.cell_name_lo; {sprintf('%d_%d',rez_mec.cell_id_lo(i,1),rez_mec.cell_id_lo(i,2))}];
end
% find intersection
rez_mec.uniq_name_hi = unique(rez_mec.cell_name_hi);
rez_mec.uniq_name_lo = unique(rez_mec.cell_name_lo);
rez_mec.uniq_name_both = intersect(rez_mec.uniq_name_hi,rez_mec.uniq_name_lo);
% take id info for cells with at least one stable rep in both high and low
rez_mec.cell_id_hi_intersect = rez_mec.cell_id_hi(ismember(rez_mec.cell_name_hi,rez_mec.uniq_name_both),:);
rez_mec.cell_id_lo_intersect = rez_mec.cell_id_lo(ismember(rez_mec.cell_name_lo,rez_mec.uniq_name_both),:);

% VISp
% dim1=sesh num, 2=cell id, 3=rep num, 4=gc trial num
rez_vis.num_cells_total = 0;
rez_vis.cell_id_hi = [];
rez_vis.cell_id_lo = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_vis.repinfo_hi)    
    if ~isempty(rez_vis.corrmat_hi{i})

        cell_id_this = rez_vis.repinfo_hi{i}.CellID{1};
        rez_vis.num_cells_total = rez_vis.num_cells_total + numel(cell_id_this);
        cell_id_hi_this = [];
        cell_id_lo_this = [];
        
        % high contrast        
        trial_corr_mat_this = rez_vis.corrmat_hi{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_vis.repinfo_hi{i},1)
                stab = rez_vis.repinfo_hi{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_hi_this = [cell_id_hi_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end       

        % low contrast
        trial_corr_mat_this = rez_vis.corrmat_lo{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_vis.repinfo_lo{i},1)
                stab = rez_vis.repinfo_lo{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_lo_this = [cell_id_lo_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end

        % add to list
        rez_vis.cell_id_hi = [rez_vis.cell_id_hi; cell_id_hi_this];
        rez_vis.cell_id_lo = [rez_vis.cell_id_lo; cell_id_lo_this];
    end
end
% create id strings for cells
rez_vis.cell_name_hi = [];
rez_vis.cell_name_lo = [];
for i = 1:size(rez_vis.cell_id_hi,1)
    rez_vis.cell_name_hi = [rez_vis.cell_name_hi; {sprintf('%d_%d',rez_vis.cell_id_hi(i,1),rez_vis.cell_id_hi(i,2))}];
end
for i = 1:size(rez_vis.cell_id_lo,1)
    rez_vis.cell_name_lo = [rez_vis.cell_name_lo; {sprintf('%d_%d',rez_vis.cell_id_lo(i,1),rez_vis.cell_id_lo(i,2))}];
end
% find intersection
rez_vis.uniq_name_hi = unique(rez_vis.cell_name_hi);
rez_vis.uniq_name_lo = unique(rez_vis.cell_name_lo);
rez_vis.uniq_name_both = intersect(rez_vis.uniq_name_hi,rez_vis.uniq_name_lo);
% take id info for cells with at least one stable rep in both high and low
rez_vis.cell_id_hi_intersect = rez_vis.cell_id_hi(ismember(rez_vis.cell_name_hi,rez_vis.uniq_name_both),:);
rez_vis.cell_id_lo_intersect = rez_vis.cell_id_lo(ismember(rez_vis.cell_name_lo,rez_vis.uniq_name_both),:);

% RSP
% dim1=sesh num, 2=cell id, 3=rep num, 4=gc trial num
rez_rsp.num_cells_total = 0;
rez_rsp.cell_id_hi = [];
rez_rsp.cell_id_lo = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_rsp.repinfo_hi)    
    if ~isempty(rez_rsp.corrmat_hi{i})

        cell_id_this = rez_rsp.repinfo_hi{i}.CellID{1};
        rez_rsp.num_cells_total = rez_rsp.num_cells_total + numel(cell_id_this);
        cell_id_hi_this = [];
        cell_id_lo_this = [];
        
        % high contrast        
        trial_corr_mat_this = rez_rsp.corrmat_hi{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_rsp.repinfo_hi{i},1)
                stab = rez_rsp.repinfo_hi{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_hi_this = [cell_id_hi_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end       

        % low contrast
        trial_corr_mat_this = rez_rsp.corrmat_lo{i};
        for ii = 1:numel(cell_id_this)
            for j = 1:size(rez_rsp.repinfo_lo{i},1)
                stab = rez_rsp.repinfo_lo{i}.Stab{j}(ii) > opt.stab_thresh;
                for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                    non_remap = nanmean(squeeze(trial_corr_mat_this(j,ii,k,1:opt.num_tr_bl))) > opt.remap_thresh;
                    if stab && non_remap
                        cell_id_lo_this = [cell_id_lo_this; [i cell_id_this(ii) j k-opt.num_tr_bl]];
                    end
                end
            end
        end

        % add to list
        rez_rsp.cell_id_hi = [rez_rsp.cell_id_hi; cell_id_hi_this];
        rez_rsp.cell_id_lo = [rez_rsp.cell_id_lo; cell_id_lo_this];
    end
end
% create id strings for cells
rez_rsp.cell_name_hi = [];
rez_rsp.cell_name_lo = [];
for i = 1:size(rez_rsp.cell_id_hi,1)
    rez_rsp.cell_name_hi = [rez_rsp.cell_name_hi; {sprintf('%d_%d',rez_rsp.cell_id_hi(i,1),rez_rsp.cell_id_hi(i,2))}];
end
for i = 1:size(rez_rsp.cell_id_lo,1)
    rez_rsp.cell_name_lo = [rez_rsp.cell_name_lo; {sprintf('%d_%d',rez_rsp.cell_id_lo(i,1),rez_rsp.cell_id_lo(i,2))}];
end
% find intersection
rez_rsp.uniq_name_hi = unique(rez_rsp.cell_name_hi);
rez_rsp.uniq_name_lo = unique(rez_rsp.cell_name_lo);
rez_rsp.uniq_name_both = intersect(rez_rsp.uniq_name_hi,rez_rsp.uniq_name_lo);
% take id info for cells with at least one stable rep in both high and low
rez_rsp.cell_id_hi_intersect = rez_rsp.cell_id_hi(ismember(rez_rsp.cell_name_hi,rez_rsp.uniq_name_both),:);
rez_rsp.cell_id_lo_intersect = rez_rsp.cell_id_lo(ismember(rez_rsp.cell_name_lo,rez_rsp.uniq_name_both),:);

%% compute firing rates around peak
opt.num_bins = 2*opt.max_lag_peak_fr/opt.SpatialBin+1;

% MEC
fr_hi = nan(numel(rez_mec.uniq_name_both),opt.num_bins);
fr_lo = nan(numel(rez_mec.uniq_name_both),opt.num_bins);
tic
parfor i = 1:numel(rez_mec.uniq_name_both)
    fprintf('computing peak aligned fr: cell %d/%d\n',i,numel(rez_mec.uniq_name_both));
    cell_id_hi_this = rez_mec.cell_id_hi(strcmp(rez_mec.cell_name_hi,rez_mec.uniq_name_both(i)),:);
    cell_id_lo_this = rez_mec.cell_id_lo(strcmp(rez_mec.cell_name_lo,rez_mec.uniq_name_both(i)),:);   
    session_name_this = rez_mec.session_name{cell_id_hi_this(1,1)};
    cell_id_this = cell_id_hi_this(1,2);
    reps_trials_hi = cell_id_hi_this(:,3:4);
    reps_trials_lo = cell_id_lo_this(:,3:4);
    dat = load(fullfile(paths.data,session_name_this));
    [fr_hi(i,:),fr_lo(i,:)] = compute_peak_aligned_fr_gaincontrast(cell_id_this,reps_trials_hi,reps_trials_lo,dat,opt);
end
toc
rez_mec.fr_hi = fr_hi;
rez_mec.fr_lo = fr_lo;

% VISp
fr_hi = nan(numel(rez_vis.uniq_name_both),opt.num_bins);
fr_lo = nan(numel(rez_vis.uniq_name_both),opt.num_bins);
tic
parfor i = 1:numel(rez_vis.uniq_name_both)
    fprintf('computing peak aligned fr: cell %d/%d\n',i,numel(rez_vis.uniq_name_both));
    cell_id_hi_this = rez_vis.cell_id_hi(strcmp(rez_vis.cell_name_hi,rez_vis.uniq_name_both(i)),:);
    cell_id_lo_this = rez_vis.cell_id_lo(strcmp(rez_vis.cell_name_lo,rez_vis.uniq_name_both(i)),:);   
    session_name_this = rez_vis.session_name{cell_id_hi_this(1,1)};
    cell_id_this = cell_id_hi_this(1,2);
    reps_trials_hi = cell_id_hi_this(:,3:4);
    reps_trials_lo = cell_id_lo_this(:,3:4);
    dat = load(fullfile(paths.data,session_name_this));
    [fr_hi(i,:),fr_lo(i,:)] = compute_peak_aligned_fr_gaincontrast(cell_id_this,reps_trials_hi,reps_trials_lo,dat,opt);
end
toc
rez_vis.fr_hi = fr_hi;
rez_vis.fr_lo = fr_lo;

% RSP
fr_hi = nan(numel(rez_rsp.uniq_name_both),opt.num_bins);
fr_lo = nan(numel(rez_rsp.uniq_name_both),opt.num_bins);
tic
parfor i = 1:numel(rez_rsp.uniq_name_both)
    fprintf('computing peak aligned fr: cell %d/%d\n',i,numel(rez_rsp.uniq_name_both));
    cell_id_hi_this = rez_rsp.cell_id_hi(strcmp(rez_rsp.cell_name_hi,rez_rsp.uniq_name_both(i)),:);
    cell_id_lo_this = rez_rsp.cell_id_lo(strcmp(rez_rsp.cell_name_lo,rez_rsp.uniq_name_both(i)),:);   
    session_name_this = rez_rsp.session_name{cell_id_hi_this(1,1)};
    cell_id_this = cell_id_hi_this(1,2);
    reps_trials_hi = cell_id_hi_this(:,3:4);
    reps_trials_lo = cell_id_lo_this(:,3:4);
    dat = load(fullfile(paths.data,session_name_this));
    [fr_hi(i,:),fr_lo(i,:)] = compute_peak_aligned_fr_gaincontrast(cell_id_this,reps_trials_hi,reps_trials_lo,dat,opt);
end
toc
rez_rsp.fr_hi = fr_hi;
rez_rsp.fr_lo = fr_lo;

%% Figure 5F

xplot = -opt.max_lag_peak_fr:opt.SpatialBin:opt.max_lag_peak_fr;
hfig = figure('Position',[200 200 1000 300]); 
hfig.Name = 'Figure 5F: peak aligned firing rate MEC VISp RSP';

% MEC
subplot(1,3,1); hold on;
means = [nanmean(rez_mec.fr_hi); nanmean(rez_mec.fr_lo)];
sems = [nanstd(rez_mec.fr_hi); nanstd(rez_mec.fr_lo)]/sqrt(size(rez_mec.fr_hi,1));
shadedErrorBar(xplot,means(1,:),sems(1,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_hi)});
shadedErrorBar(xplot,means(2,:),sems(2,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_lo)});
xlim([-opt.max_lag_peak_fr opt.max_lag_peak_fr]);
%ylim([0.1 0.9]);
xticks(-opt.max_lag_peak_fr:10:opt.max_lag_peak_fr);
%yticks(0.1:0.2:0.9);
plot([0 0],ylim(),'k--');
ylabel('Firing rate');
title(sprintf('MEC (%d units)',size(rez_mec.fr_hi,1)));
legend(sprintf('c=%d',opt.contr_hi),sprintf('c=%d',opt.contr_lo))

% VISp
subplot(1,3,2); hold on;
means = [nanmean(rez_vis.fr_hi); nanmean(rez_vis.fr_lo)];
sems = [nanstd(rez_vis.fr_hi); nanstd(rez_vis.fr_lo)]/sqrt(size(rez_vis.fr_hi,1));
shadedErrorBar(xplot,means(1,:),sems(1,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_hi)});
shadedErrorBar(xplot,means(2,:),sems(2,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_lo)});
xlabel('Dist. from baseline peak (cm)');
title(sprintf('VISp (%d units)',size(rez_vis.fr_hi,1)));
xlim([-opt.max_lag_peak_fr opt.max_lag_peak_fr]);
%ylim([0.1 0.9]);
xticks(-opt.max_lag_peak_fr:10:opt.max_lag_peak_fr);
yticks(0.1:0.2:0.9);
plot([0 0],ylim(),'k--');

% plot: RSP
subplot(1,3,3); hold on;
means = [nanmean(rez_rsp.fr_hi); nanmean(rez_rsp.fr_lo)];
sems = [nanstd(rez_rsp.fr_hi); nanstd(rez_rsp.fr_lo)]/sqrt(size(rez_rsp.fr_hi,1));
shadedErrorBar(xplot,means(1,:),sems(1,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_hi)});
shadedErrorBar(xplot,means(2,:),sems(2,:),'lineProps',{'Color',get_color(opt.gain,opt.contr_lo)});
title(sprintf('RSP (%d units)',size(rez_rsp.fr_hi,1)));
xlim([-opt.max_lag_peak_fr opt.max_lag_peak_fr]);
%ylim([0.1 0.9]);
xticks(-opt.max_lag_peak_fr:10:opt.max_lag_peak_fr);
%yticks(0.1:0.2:0.9);
plot([0 0],ylim(),'k--');
