% Makes the following figures:
% Fig 6P
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.sessions = fullfile(paths.intermediate_data,'session_lists');
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

% opts for gain contrast analysis
opt.brain_region = 'MEC';
opt.dark = false;
opt.max_lag = 30;
opt.gain = 0.8;
opt.contr_hi = 100;
opt.contr_lo_all = [10 20];
opt.num_tr_bl = 6;
opt.stab_thresh = 0.5;
opt.remap_thresh = 0.5;
opt.min_num_stab_cells = 5; 

%% iterate over contrasts
y_all = {};
grp = {};
for ctrIdx = 1:numel(opt.contr_lo_all)
    
    opt.contr_lo = opt.contr_lo_all(ctrIdx);
    
    session_file = sprintf('small_gain_contrast%d_MEC',opt.contr_lo);
    rez_mec = struct; % struct for holding MEC results
    [rez_mec.corrmat_hi, rez_mec.corrmat_lo, rez_mec.shiftmat_hi, rez_mec.shiftmat_lo, rez_mec.repinfo_hi, rez_mec.repinfo_lo] = ...
        xcorr_peak_and_shift_gaincontrast(session_file,cell_info,opt,paths);

    rez_mec.all_shifts = [];
    rez_mec.sesh_nums = [];
    rez_mec.cell_ids = [];
    rez_mec.dist_tuned = [];
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



            % get cell ids
            cell_id_this = rez_mec.repinfo_hi{i}.CellID{1};
            session_this = rez_mec.repinfo_hi{i}.Session{1};
            session_split = strsplit(session_this,'_');
            cell_id_full = cell(numel(cell_id_this),1);
            dist_tuned_this = nan(numel(cell_id_this),1);
            for j = 1:numel(cell_id_this)
                cell_id_full{j} = sprintf('%s_%s_c%d',session_split{1},session_split{2},cell_id_this(j));
                tmp = dist_tuned(strcmp(dt.UniqueID,cell_id_full{j}));
                if ~isempty(tmp)
                    dist_tuned_this(j) = tmp;
                end
            end

            % both
            keep = ~isnan(avg_shift_by_cell_hi) & ~isnan(avg_shift_by_cell_lo);
            rez_mec.all_shifts = [rez_mec.all_shifts; [avg_shift_by_cell_hi(keep) avg_shift_by_cell_lo(keep)]];
            rez_mec.sesh_nums = [rez_mec.sesh_nums; i*ones(sum(keep),1)];
            rez_mec.cell_ids = [rez_mec.cell_ids; cell_id_full(keep)];
            rez_mec.dist_tuned = [rez_mec.dist_tuned; dist_tuned_this(keep)];

        end
    end
    
    y_all{1+(ctrIdx-1)*2} = y1;
    y_all{2+(ctrIdx-1)*2} = y2;
end

%% Figure 6P: ANOVA on contrast, dist vs non-dist

y_concat = [y_all{1}; y_all{3}; y_all{2}; y_all{4}];
grp = {};
grp{1} = [zeros(numel(y_all{1})+numel(y_all{3}),1); ones(numel(y_all{2})+numel(y_all{4}),1)];
grp{2} = [zeros(size(y_all{1})); ones(size(y_all{3})); zeros(size(y_all{2})); ones(size(y_all{4}))];

[anovap,t,stats] = anovan(y_concat,grp);
grp_all = [ones(size(y_all{1})); 2*ones(size(y_all{3})); 3*ones(size(y_all{2})); 4*ones(size(y_all{4}))];

hfig = figure('Position',[200 200 250 450]); hold on;
hfig.Name = 'Figure 6P: change in map shift_dist vs non dist c10 and c20';
boxplot(y_concat,grp_all);
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','o'); % Change symbols for all the groups.
ylabel('Change in shift (cm)');
plot(xlim(),[0 0],'k--');
xticklabels({'Dist, c=10','Dist, c=20','Non-Dist, c=10','Non-Dist, c=20'});
xtickangle(90);
title(sprintf('ANOVA p, Dist vs Non-Dist = %0.3f',anovap(1)));