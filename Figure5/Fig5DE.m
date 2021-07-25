% Makes the following figures:
% Fig 5D,E
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
opt.min_num_stab_cells = 5; 

session_file = {'small_gain_contrast10_MEC','small_gain_contrast10_VISp','small_gain_contrast10_RSP'};

%% MEC
opt.brain_region = 'MEC';
rez_mec = struct; % struct for holding MEC results
[rez_mec.corrmat_hi, rez_mec.corrmat_lo, rez_mec.shiftmat_hi, rez_mec.shiftmat_lo, rez_mec.repinfo_hi, rez_mec.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(session_file{1},cell_info,opt,paths);

% VISp
opt.brain_region = 'VISp';
rez_vis = struct; % struct for holding vis results
[rez_vis.corrmat_hi, rez_vis.corrmat_lo, rez_vis.shiftmat_hi, rez_vis.shiftmat_lo, rez_vis.repinfo_hi, rez_vis.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(session_file{2},cell_info,opt,paths);

% RSP
opt.brain_region = 'RSP';
rez_rsp = struct; % struct for holding vis results
[rez_rsp.corrmat_hi, rez_rsp.corrmat_lo, rez_rsp.shiftmat_hi, rez_rsp.shiftmat_lo, rez_rsp.repinfo_hi, rez_rsp.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(session_file{3},cell_info,opt,paths);

%% Extract map shifts for high and low contrast

% MEC
rez_mec.all_shifts = [];
rez_mec.sesh_nums = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_mec.repinfo_hi)    
    if ~isempty(rez_mec.corrmat_hi{i})
        % high contrast
        trial_corr_mat_this = rez_mec.corrmat_hi{i};
        shift_mat_this = rez_mec.shiftmat_hi{i};
        for j = 1:size(rez_mec.repinfo_hi{i},1)
            stab_cells = rez_mec.repinfo_hi{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_hi = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_hi = squeeze(nanmean(avg_shift_by_cell_hi))';

        % low contrast
        trial_corr_mat_this = rez_mec.corrmat_lo{i};
        shift_mat_this = rez_mec.shiftmat_lo{i};
        for j = 1:size(rez_mec.repinfo_lo{i},1)
            stab_cells = rez_mec.repinfo_lo{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_lo = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_lo = squeeze(nanmean(avg_shift_by_cell_lo))';

        % both
        keep = ~isnan(avg_shift_by_cell_hi) & ~isnan(avg_shift_by_cell_lo);
        rez_mec.all_shifts = [rez_mec.all_shifts; [avg_shift_by_cell_hi(keep) avg_shift_by_cell_lo(keep)]];
        rez_mec.sesh_nums = [rez_mec.sesh_nums; i*ones(sum(keep),1)];
    end
end

% VISp
rez_vis.all_shifts = [];
rez_vis.sesh_nums = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_vis.repinfo_hi)    
    if ~isempty(rez_vis.corrmat_hi{i})
        % high contrast
        trial_corr_mat_this = rez_vis.corrmat_hi{i};
        shift_mat_this = rez_vis.shiftmat_hi{i};
        for j = 1:size(rez_vis.repinfo_hi{i},1)
            stab_cells = rez_vis.repinfo_hi{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_hi = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_hi = squeeze(nanmean(avg_shift_by_cell_hi))';

        % low contrast
        trial_corr_mat_this = rez_vis.corrmat_lo{i};
        shift_mat_this = rez_vis.shiftmat_lo{i};
        for j = 1:size(rez_vis.repinfo_lo{i},1)
            stab_cells = rez_vis.repinfo_lo{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_lo = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_lo = squeeze(nanmean(avg_shift_by_cell_lo))';

        % both
        keep = ~isnan(avg_shift_by_cell_hi) & ~isnan(avg_shift_by_cell_lo);
        rez_vis.all_shifts = [rez_vis.all_shifts; [avg_shift_by_cell_hi(keep) avg_shift_by_cell_lo(keep)]];
        rez_vis.sesh_nums = [rez_vis.sesh_nums; i*ones(sum(keep),1)];
    end
end

% RSP
rez_rsp.all_shifts = [];
rez_rsp.sesh_nums = [];
% find cells with at least one stable trial in both hi and lo contrast
for i = 1:numel(rez_rsp.repinfo_hi)    
    if ~isempty(rez_rsp.corrmat_hi{i})
        % high contrast
        trial_corr_mat_this = rez_rsp.corrmat_hi{i};
        shift_mat_this = rez_rsp.shiftmat_hi{i};
        for j = 1:size(rez_rsp.repinfo_hi{i},1)
            stab_cells = rez_rsp.repinfo_hi{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_hi = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_hi = squeeze(nanmean(avg_shift_by_cell_hi))';

        % low contrast
        trial_corr_mat_this = rez_rsp.corrmat_lo{i};
        shift_mat_this = rez_rsp.shiftmat_lo{i};
        for j = 1:size(rez_rsp.repinfo_lo{i},1)
            stab_cells = rez_rsp.repinfo_lo{i}.Stab{j} > opt.stab_thresh;
            shift_mat_this(j,~stab_cells,:,:) = nan;
            for k = opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc
                non_remap_cells = nanmean(squeeze(trial_corr_mat_this(j,:,k,1:opt.num_tr_bl)),2) > opt.remap_thresh;
                shift_mat_this(j,~non_remap_cells,k,:) = nan;
                shift_mat_this(j,~non_remap_cells,:,k) = nan;
            end
        end
        avg_shift_by_cell_lo = squeeze(nanmean(nanmean(shift_mat_this(:,:,1:opt.num_tr_bl,opt.num_tr_bl+1:opt.num_tr_bl+opt.num_tr_gc),4),3));
        avg_shift_by_cell_lo = squeeze(nanmean(avg_shift_by_cell_lo))';

        % both
        keep = ~isnan(avg_shift_by_cell_hi) & ~isnan(avg_shift_by_cell_lo);
        rez_rsp.all_shifts = [rez_rsp.all_shifts; [avg_shift_by_cell_hi(keep) avg_shift_by_cell_lo(keep)]];
        rez_rsp.sesh_nums = [rez_rsp.sesh_nums; i*ones(sum(keep),1)];
    end
end

%% Figure 5D: violin plot of shift change
shift_change_mec = rez_mec.all_shifts(:,2)-rez_mec.all_shifts(:,1);
shift_change_vis = rez_vis.all_shifts(:,2)-rez_vis.all_shifts(:,1);
shift_change_rsp = rez_rsp.all_shifts(:,2)-rez_rsp.all_shifts(:,1);

hfig = figure('Position',[400 400 200 300]); hold on; 
hfig.Name = 'Figure 5D: change in map shift_brain region comparison';
distributionPlot({shift_change_mec,shift_change_vis,shift_change_rsp},'color',{[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5]},'showMM',6);
% boxplot([shift_change_mec; shift_change_vis; shift_change_rsp],...
%     [ones(numel(shift_change_mec),1); 2*ones(numel(shift_change_vis),1); 3*ones(numel(shift_change_rsp),1)]);
xticklabels({'MEC','V1','RSC'});
ylabel('Change in shift (cm)');
plot(xlim(),[0 0],'k--');
ylim([-30 30]);
% statistical test

pvals = [];
pvals(1) = signrank(shift_change_mec);
pvals(2) = signrank(shift_change_vis);
pvals(3) = signrank(shift_change_rsp);
pvals(4) = ranksum(shift_change_mec,shift_change_vis);
pvals(5) = ranksum(shift_change_mec,shift_change_rsp);
pvals(6) = ranksum(shift_change_vis,shift_change_rsp);
for i = 1:3
    if pvals(i) < 0.001
        text(i,max(ylim),'***','HorizontalAlignment','center');
    elseif pvals(i) < 0.01
        text(i,max(ylim),'**','HorizontalAlignment','center');
    elseif pvals(i) < 0.05
        text(i,max(ylim),'*','HorizontalAlignment','center');
    end
end


%% get median shifts by session
rez_mec.median_shift_diff = nan(1,numel(rez_mec.repinfo_hi));
rez_mec.median_shift_diff_pval = nan(1,numel(rez_mec.repinfo_hi));
for i = 1:numel(rez_mec.repinfo_hi)
    if sum(rez_mec.sesh_nums==i)>=opt.min_num_stab_cells
         shift_diff_this = rez_mec.all_shifts(rez_mec.sesh_nums==i,2)-rez_mec.all_shifts(rez_mec.sesh_nums==i,1);
         rez_mec.median_shift_diff(i) = median(shift_diff_this);
         rez_mec.median_shift_diff_pval(i) = signrank(shift_diff_this);
    end
end
keep = ~isnan(rez_mec.median_shift_diff);
rez_mec.median_shift_diff = rez_mec.median_shift_diff(keep);
rez_mec.median_shift_diff_pval = rez_mec.median_shift_diff_pval(keep);

rez_vis.median_shift_diff = nan(1,numel(rez_vis.repinfo_hi));
rez_vis.median_shift_diff_pval = nan(1,numel(rez_vis.repinfo_hi));
for i = 1:numel(rez_vis.repinfo_hi)
    if sum(rez_vis.sesh_nums==i)>=opt.min_num_stab_cells
         shift_diff_this = rez_vis.all_shifts(rez_vis.sesh_nums==i,2)-rez_vis.all_shifts(rez_vis.sesh_nums==i,1);
         rez_vis.median_shift_diff(i) = median(shift_diff_this);
         rez_vis.median_shift_diff_pval(i) = signrank(shift_diff_this);
    end
end
keep = ~isnan(rez_vis.median_shift_diff);
rez_vis.median_shift_diff = rez_vis.median_shift_diff(keep);
rez_vis.median_shift_diff_pval = rez_vis.median_shift_diff_pval(keep);

rez_rsp.median_shift_diff = nan(1,numel(rez_rsp.repinfo_hi));
rez_rsp.median_shift_diff_pval = nan(1,numel(rez_rsp.repinfo_hi));
for i = 1:numel(rez_rsp.repinfo_hi)
    if sum(rez_rsp.sesh_nums==i)>=opt.min_num_stab_cells
         shift_diff_this = rez_rsp.all_shifts(rez_rsp.sesh_nums==i,2)-rez_rsp.all_shifts(rez_rsp.sesh_nums==i,1);
         rez_rsp.median_shift_diff(i) = median(shift_diff_this);
         rez_rsp.median_shift_diff_pval(i) = signrank(shift_diff_this);
    end
end
keep = ~isnan(rez_rsp.median_shift_diff);
rez_rsp.median_shift_diff = rez_rsp.median_shift_diff(keep);
rez_rsp.median_shift_diff_pval = rez_rsp.median_shift_diff_pval(keep);

%% Figure 5E: plot session medians with size + fill color determined by p value

hfig = figure('Position',[400 400 200 300]); hold on; 
hfig.Name = 'Figure 5E: change in map shift_session median_brain region comparison';

% MEC
for i = 1:numel(rez_mec.median_shift_diff)
    jitt = 0.03*randn(1,1);
    if rez_mec.median_shift_diff_pval(i)<0.05
        myscatter(1+jitt,rez_mec.median_shift_diff(i),'r',0.5);
    else
        myscatter(1+jitt,rez_mec.median_shift_diff(i),'k',0.2);
    end
end

% VISp
for i = 1:numel(rez_vis.median_shift_diff)
    jitt = 0.03*randn(1,1);
    if rez_vis.median_shift_diff_pval(i)<0.05
        myscatter(2+jitt,rez_vis.median_shift_diff(i),'r',0.5);
    else
        myscatter(2+jitt,rez_vis.median_shift_diff(i),'k',0.2);
    end
end

% RSP
for i = 1:numel(rez_rsp.median_shift_diff)
    jitt = 0.03*randn(1,1);
    if rez_rsp.median_shift_diff_pval(i)<0.05
        myscatter(3+jitt,rez_rsp.median_shift_diff(i),'r',0.5);
    else
        myscatter(3+jitt,rez_rsp.median_shift_diff(i),'k',0.2);
    end
end

xticks([1 2 3]);
xticklabels({'MEC','V1','RSC'})
xlim([0.5 3.5]);
ylim([-8 6]);
ylabel('Change in shift (cm)');
plot(xlim(),[0 0],'k:');

pval_median(1) = signrank(rez_mec.median_shift_diff);
pval_median(2) = signrank(rez_vis.median_shift_diff);
pval_median(3) = signrank(rez_rsp.median_shift_diff);
for i = 1:3
    if pval_median(i) < 0.001
        text(i,max(ylim),'***','HorizontalAlignment','center');
    elseif pval_median(i) < 0.01
        text(i,max(ylim),'**','HorizontalAlignment','center');
    elseif pval_median(i) < 0.05
        text(i,max(ylim),'*','HorizontalAlignment','center');
    end
end

