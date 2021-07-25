% Makes the following figures:
% Fig S4H
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% MGC 7/19/2021

%restoredefaultpath;

% add helper functions to path
%addpath(genpath('matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.sessions = fullfile(paths.intermediate_data,'session_lists');

% load cell_info table
load(fullfile(paths.intermediate_data,fullfile('cell_info','cell_info_Apr2020'))); 

% analysis options
opt = load_default_opt;
opt.gain = 0.8;
opt.contr_hi = 100;
opt.contr_lo = 20;
opt.num_tr_bl = 6;
opt.stab_thresh = 0.5;
opt.remap_thresh = 0.5;
opt.min_num_stab_cells = 5; 

%% MEC
opt.brain_region = 'MEC';
session_file = 'small_gain_contrast20_MEC';
rez_mec = struct; % struct for holding MEC results
[rez_mec.corrmat_hi, rez_mec.corrmat_lo, rez_mec.shiftmat_hi, rez_mec.shiftmat_lo, rez_mec.repinfo_hi, rez_mec.repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(session_file,cell_info,opt,paths);

%% MEC
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

%% Figure S4H: Violin plot of shift change
shift_change_mec = rez_mec.all_shifts(:,2)-rez_mec.all_shifts(:,1);
hfig = figure('Position',[400 400 350 400]); hold on; 
hfig.Name = 'Figure S4H: change in map shift_brain region comparison';
distributionPlot({shift_change_mec},'color',{[0.5 0.5 0.5]},'showMM',6);
xticklabels({sprintf('MEC (n=%d)',size(shift_change_mec,1))});
ylabel('change in shift (cm)');
plot(xlim(),[0 0],'k--');
ylim([-25 25]);

% statistical test
pval = signrank(shift_change_mec);