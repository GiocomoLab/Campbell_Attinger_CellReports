function [corrmat,shiftmat,rep_info,corrmat_bl,shiftmat_bl,rep_info_bl] = compute_trial_corr_mat_indiv_neurons_with_baseline(session_name,cell_info,paths,opt)
% compute_trial_corr_mat_indiv_neurons.m 
%   MGC 3/30/2020
% remember that this takes max corr over -opt.max_lag:opt.max_lag

corrmat_all = {};
shiftmat_all = {};
rep_info_all = {};
corrmat_all_bl = {};
shiftmat_all_bl = {};
rep_info_all_bl = {};
pb = ParforProgressbar(numel(session_name));
tic
parfor i = 1:numel(session_name)
    % fprintf('analyzing %s session %d/%d: %s\n',opt.brain_region,i,numel(session_name),session_name{i});
    if ~exist(fullfile(paths.data,sprintf('%s.mat',session_name{i})),'file')
        fprintf('WARNING: %s not found\n',session_name{i})
        continue
    end
    dat = load(fullfile(paths.data,session_name{i}));
    if strcmp(opt.brain_region,'RSP')
        good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{i}) & ...
            (strcmp(cell_info.BrainRegion,'RSPd') | strcmp(cell_info.BrainRegion,'RSPv'))); 
    else
        good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{i}) & ...
            strcmp(cell_info.BrainRegion,opt.brain_region));
    end         
    mouse = session_name{i}(1:4);
    session = session_name(i);
    
    % find blocks of high contrast gain trials
    trials = find(dat.trial_gain==opt.gain & dat.trial_contrast==100);
    trials = reshape(trials,4,numel(trials)/4)';
    trials2 = nan(size(trials,1),size(trials,2)+2*opt.num_tr_bl);
    for j = 1:size(trials,1)
        trials2(j,:) = trials(j,1)-opt.num_tr_bl:trials(j,end)+opt.num_tr_bl;
    end
    trials = trials2(sum(dat.trial_contrast(trials2')~=100)==0,:); % make sure all high contrast
    
    % cell info for this session
    N = size(trials,1)*numel(good_cells);
    cell_info_this = table;
    cell_info_this.Mouse = repmat(mouse,N,1);
    cell_info_this.Session = repmat(session,N,1);
    cell_info_this.CellID = repmat(good_cells,size(trials,1),1);
    cell_info_this.Rep = reshape(repmat(1:size(trials,1),numel(good_cells),1),N,1);
    cell_info_this.Trials = nan(N,size(trials,2));
    
    % compute trial corr mat for each rep, keeping only stable cells
    trial_corr_mat_this = nan(N,size(trials,2),size(trials,2));
    shift_mat_this = nan(N,size(trials,2),size(trials,2));
    for j = 1:size(trials,1)
        [xcorr_tmp,~,shift_tmp] = trialCorrMat(good_cells,trials(j,:),dat,opt);
        trial_corr_mat_this(numel(good_cells)*(j-1)+1:numel(good_cells)*j,:,:) = xcorr_tmp;
        shift_mat_this(numel(good_cells)*(j-1)+1:numel(good_cells)*j,:,:) = shift_tmp;
        cell_info_this.Trials(numel(good_cells)*(j-1)+1:numel(good_cells)*j,:) = repmat(trials(j,:),numel(good_cells),1);
    end
    stab_this = nanmean(nanmean(trial_corr_mat_this(:,1:opt.num_tr_bl,1:opt.num_tr_bl),3),2);
    cell_info_this.Stab_BL = stab_this;
    cell_info_this.Gain = dat.trial_gain(cell_info_this.Trials);
    
    % save data for later
    corrmat_all{i} = trial_corr_mat_this;
    shiftmat_all{i} = shift_mat_this;
    rep_info_all{i} = cell_info_this;
    
    %% REPEAT FOR BASELINE TRIALS
    
    % find preceding blocks of high contrast baseline trials
    trials = trials-opt.num_tr_bl-opt.num_tr_gc;
    trials = trials(sum(dat.trial_contrast(trials')~=100)==0 & ...
        sum(dat.trial_gain(trials')~=1)==0,:); % make sure all high contrast and gain=1
    
    % cell info for this session
    N = size(trials,1)*numel(good_cells);
    cell_info_this = table;
    cell_info_this.Mouse = repmat(mouse,N,1);
    cell_info_this.Session = repmat(session,N,1);
    cell_info_this.CellID = repmat(good_cells,size(trials,1),1);
    cell_info_this.Rep = reshape(repmat(1:size(trials,1),numel(good_cells),1),N,1);
    
    % compute trial corr mat for each rep, keeping only stable cells
    trial_corr_mat_this = nan(N,size(trials,2),size(trials,2));
    shift_mat_this = nan(N,size(trials,2),size(trials,2));
    for j = 1:size(trials,1)
        [xcorr_tmp,~,shift_tmp] = trialCorrMat(good_cells,trials(j,:),dat,opt);
        trial_corr_mat_this(numel(good_cells)*(j-1)+1:numel(good_cells)*j,:,:) = xcorr_tmp;
        shift_mat_this(numel(good_cells)*(j-1)+1:numel(good_cells)*j,:,:) = shift_tmp;
    end
    stab_this = nanmean(nanmean(trial_corr_mat_this(:,1:opt.num_tr_bl,1:opt.num_tr_bl),3),2);
    cell_info_this.Stab_BL = stab_this;
    
    % save data for later
    corrmat_all_bl{i} = trial_corr_mat_this;
    shiftmat_all_bl{i} = shift_mat_this;
    rep_info_all_bl{i} = cell_info_this;
    
    pb.increment();
end
toc;

%% concat all cells

% gain change
corrmat = [];
shiftmat = [];
rep_info = [];
for i = 1:numel(corrmat_all)
    corrmat = [corrmat; corrmat_all{i}];
    shiftmat = [shiftmat; shiftmat_all{i}];
    rep_info = [rep_info; rep_info_all{i}];
end
rep_info.rep_id = strcat(rep_info.Session,'_rep',num2str(rep_info.Rep));

% baseline
corrmat_bl = [];
shiftmat_bl = [];
rep_info_bl = [];
for i = 1:numel(corrmat_all_bl)
    corrmat_bl = [corrmat_bl; corrmat_all_bl{i}];
    shiftmat_bl = [shiftmat_bl; shiftmat_all_bl{i}];
    rep_info_bl = [rep_info_bl; rep_info_all_bl{i}];
end
rep_info_bl.rep_id = strcat(rep_info_bl.Session,'_rep',num2str(rep_info_bl.Rep));

end

