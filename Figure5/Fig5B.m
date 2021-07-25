% Makes the following figures:
% Fig 5B
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

%% proportion of cells

% MEC
rez_mec.cell_counts = nan(1,4);
rez_mec.cell_counts(1) = numel(rez_mec.uniq_name_hi)-numel(rez_mec.uniq_name_both);
rez_mec.cell_counts(2) = numel(rez_mec.uniq_name_both);
rez_mec.cell_counts(3) = numel(rez_mec.uniq_name_lo)-numel(rez_mec.uniq_name_both);
rez_mec.cell_counts(4) = rez_mec.num_cells_total-numel(rez_mec.uniq_name_hi)-numel(rez_mec.uniq_name_lo)+numel(rez_mec.uniq_name_both);

% VIS
rez_vis.cell_counts = nan(1,4);
rez_vis.cell_counts(1) = numel(rez_vis.uniq_name_hi)-numel(rez_vis.uniq_name_both);
rez_vis.cell_counts(2) = numel(rez_vis.uniq_name_both);
rez_vis.cell_counts(3) = numel(rez_vis.uniq_name_lo)-numel(rez_vis.uniq_name_both);
rez_vis.cell_counts(4) = rez_vis.num_cells_total-numel(rez_vis.uniq_name_hi)-numel(rez_vis.uniq_name_lo)+numel(rez_vis.uniq_name_both);

% RSP
rez_rsp.cell_counts = nan(1,4);
rez_rsp.cell_counts(1) = numel(rez_rsp.uniq_name_hi)-numel(rez_rsp.uniq_name_both);
rez_rsp.cell_counts(2) = numel(rez_rsp.uniq_name_both);
rez_rsp.cell_counts(3) = numel(rez_rsp.uniq_name_lo)-numel(rez_rsp.uniq_name_both);
rez_rsp.cell_counts(4) = rez_rsp.num_cells_total-numel(rez_rsp.uniq_name_hi)-numel(rez_rsp.uniq_name_lo)+numel(rez_rsp.uniq_name_both);

%% Figure 5B
hfig = figure('Position',[100 100 600 200]);
hfig.Name = 'Figure 5B: pie chart MEC VISp RSP';
explode = [0 1 0 0];

% MEC
subplot(1,3,1);
hpie = pie(rez_mec.cell_counts, explode);
title('MEC');
hpie(1).FaceColor = 'k';
hpie(2).String = strcat(num2str(rez_mec.cell_counts(1)),' (',hpie(2).String,')');
hpie(3).FaceColor = 'g';
hpie(4).String = strcat(num2str(rez_mec.cell_counts(2)),' (',hpie(4).String,')');
hpie(5).FaceColor = get_color(1,opt.contr_lo);
hpie(6).String = strcat(num2str(rez_mec.cell_counts(3)),' (',hpie(6).String,')');
hpie(7).FaceColor = 'w';
hpie(8).String = strcat(num2str(rez_mec.cell_counts(4)),' (',hpie(8).String,')');

% VIS
subplot(1,3,2);
hpie = pie(rez_vis.cell_counts, explode);
title('VISp');
hpie(1).FaceColor = 'k';
hpie(2).String = strcat(num2str(rez_vis.cell_counts(1)),' (',hpie(2).String,')');
hpie(3).FaceColor = 'g';
hpie(4).String = strcat(num2str(rez_vis.cell_counts(2)),' (',hpie(4).String,')');
hpie(5).FaceColor = get_color(1,opt.contr_lo);
hpie(6).String = strcat(num2str(rez_vis.cell_counts(3)),' (',hpie(6).String,')');
hpie(7).FaceColor = 'w';
hpie(8).String = strcat(num2str(rez_vis.cell_counts(4)),' (',hpie(8).String,')');

% RSP
subplot(1,3,3);
hpie = pie(rez_rsp.cell_counts, explode);
title('RSP');
hpie(1).FaceColor = 'k';
hpie(2).String = strcat(num2str(rez_rsp.cell_counts(1)),' (',hpie(2).String,')');
hpie(3).FaceColor = 'g';
hpie(4).String = strcat(num2str(rez_rsp.cell_counts(2)),' (',hpie(4).String,')');
hpie(5).FaceColor = get_color(1,opt.contr_lo);
hpie(6).String = strcat(num2str(rez_rsp.cell_counts(3)),' (',hpie(6).String,')');
hpie(7).FaceColor = 'w';
hpie(8).String = strcat(num2str(rez_rsp.cell_counts(4)),' (',hpie(8).String,')');

legend({'c=100 only','both','c=10 only','none'},'Location','north');
