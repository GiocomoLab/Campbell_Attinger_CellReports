% Compute trial corr matrix (used in multiple figures)
% MGC 7/19/2021

function rez_all = compute_trial_corr_mat_all_brain_regions(paths,opt,cell_info)

rez_all = {};
for br = 1:numel(opt.brain_regions)
    rez_this = struct;
    opt.brain_region = opt.brain_regions{br};
    rez_this.brain_region = opt.brain_region;
    rez_this.session_name = read_session_file(opt.session_files{br},paths);
    [rez_this.corrmat,rez_this.shiftmat,rez_this.rep_info,...
        rez_this.corrmat_bl,rez_this.shiftmat_bl,rez_this.rep_info_bl] = ...
        compute_trial_corr_mat_indiv_neurons_with_baseline(rez_this.session_name,cell_info,paths,opt);
    
    % GAIN CHANGE
    % average by rep 
    rez_this.uniq_rep_id = unique(rez_this.rep_info.rep_id);
    rez_this.corrmat_avg = nan(numel(rez_this.uniq_rep_id),size(rez_this.corrmat,2),size(rez_this.corrmat,2));
    rez_this.shiftmat_avg = nan(numel(rez_this.uniq_rep_id),size(rez_this.corrmat,2),size(rez_this.corrmat,2));
    stab_reps = false(size(rez_this.corrmat_avg,1),1);
    for i = 1:numel(rez_this.uniq_rep_id)
        keep_cell1 = strcmp(rez_this.rep_info.rep_id,rez_this.uniq_rep_id{i});
        keep_cell2 = rez_this.rep_info.Stab_BL>opt.stab_thresh;
        keep_cell = keep_cell1 & keep_cell2;
        rez_this.corrmat_avg(i,:,:) = squeeze(nanmean(rez_this.corrmat(keep_cell,:,:),1));
        rez_this.shiftmat_avg(i,:,:) = squeeze(nanmean(rez_this.shiftmat(keep_cell,:,:),1));
        stab_reps(i) = sum(keep_cell)>=opt.min_num_stab_cells && sum(keep_cell)/sum(keep_cell1)>=opt.min_frac_stab_cells;
    end
    % only keep stable reps
    rez_this.rep_id_stab = rez_this.uniq_rep_id(stab_reps);
    rez_this.corrmat_avg_stab = rez_this.corrmat_avg(stab_reps,:,:);
    rez_this.shiftmat_avg_stab = rez_this.shiftmat_avg(stab_reps,:,:);
    
    % BASELINE
    % average by rep 
    rez_this.uniq_rep_id_bl = unique(rez_this.rep_info_bl.rep_id);
    rez_this.corrmat_avg_bl = nan(numel(rez_this.uniq_rep_id_bl),size(rez_this.corrmat_bl,2),size(rez_this.corrmat_bl,2));
    rez_this.shiftmat_avg_bl = nan(numel(rez_this.uniq_rep_id_bl),size(rez_this.corrmat_bl,2),size(rez_this.corrmat_bl,2));
    stab_reps = false(size(rez_this.corrmat_avg_bl,1),1);
    for i = 1:numel(rez_this.uniq_rep_id_bl)
        keep_cell1 = strcmp(rez_this.rep_info_bl.rep_id,rez_this.uniq_rep_id_bl{i});
        keep_cell2 = rez_this.rep_info_bl.Stab_BL>opt.stab_thresh;
        keep_cell = keep_cell1 & keep_cell2;
        rez_this.corrmat_avg_bl(i,:,:) = squeeze(nanmean(rez_this.corrmat_bl(keep_cell,:,:),1));
        rez_this.shiftmat_avg_bl(i,:,:) = squeeze(nanmean(rez_this.shiftmat_bl(keep_cell,:,:),1));
        stab_reps(i) = sum(keep_cell)>=opt.min_num_stab_cells && sum(keep_cell)/sum(keep_cell1)>=opt.min_frac_stab_cells;
    end
    % only keep stable reps
    rez_this.rep_id_stab_bl = rez_this.uniq_rep_id_bl(stab_reps);
    rez_this.corrmat_avg_stab_bl = rez_this.corrmat_avg_bl(stab_reps,:,:);
    rez_this.shiftmat_avg_stab_bl = rez_this.shiftmat_avg_bl(stab_reps,:,:);
    
    % save to big list
    rez_all{br} = rez_this;
end

end
