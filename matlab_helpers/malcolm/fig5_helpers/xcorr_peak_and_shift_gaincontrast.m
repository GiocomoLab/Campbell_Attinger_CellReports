function [corrmat_hi, corrmat_lo, shiftmat_hi, shiftmat_lo, repinfo_hi, repinfo_lo] = ...
    xcorr_peak_and_shift_gaincontrast(session_file,cell_info,opt,paths)
% makes trial-trial corr matrices for small gain sessions, hi + lo contrast
% MGC 3/12/2020

session_name = read_session_file(session_file,paths);

%% compute trial corr matrix for each session
% remember that this takes max corr over -opt.max_lag:opt.max_lag
corrmat_hi = {};
corrmat_lo = {};
shiftmat_hi = {};
shiftmat_lo = {};
repinfo_hi = {};
repinfo_lo = {};
%pb = ParforProgressbar(numel(session_name));
tic
parfor i = 1:numel(session_name)
    fprintf('analyzing session %d/%d: %s\n',i,numel(session_name),session_name{i});
    dat = load(fullfile(paths.data,session_name{i}));
    if strcmp(opt.brain_region,'RSP')
        good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{i}) & ...
            (strcmp(cell_info.BrainRegion,'RSPd') | strcmp(cell_info.BrainRegion,'RSPv')));  
    else
        good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{i}) & ...
            strcmp(cell_info.BrainRegion,opt.brain_region));  
    end
    
    %% find blocks of trials with gain change + high contrast
    trials = find(dat.trial_gain==opt.gain & dat.trial_contrast==opt.contr_hi);
    trials = reshape(trials,4,numel(trials)/4)';
    trials2 = nan(size(trials,1),size(trials,2)+2*opt.num_tr_bl);
    for j = 1:size(trials,1)
        trials2(j,:) = trials(j,1)-opt.num_tr_bl:trials(j,end)+opt.num_tr_bl;
    end
    trials = trials2(sum(dat.trial_contrast(trials2')~=opt.contr_hi)==0,:); % make sure all correct contrast
    
    % make rep info table for this session
    rep_info_this = table;
    mouse = strsplit(session_name{i},'_');
    rep_info_this.Mouse = repmat(mouse(1),size(trials,1),1);
    rep_info_this.Session = repmat(session_name(i),size(trials,1),1);
    rep_info_this.Rep = (1:size(trials,1))';
    rep_info_this.CellID = repmat({good_cells},size(trials,1),1);
    rep_info_this.Stab = repmat({nan(numel(good_cells),1)},size(trials,1),1);   
    
    % compute trial corr mat for each rep
    trial_corr_mat_this = nan(size(trials,1),numel(good_cells),size(trials,2),size(trials,2));
    shift_mat_this = nan(size(trials,1),numel(good_cells),size(trials,2),size(trials,2));
    for j = 1:size(trials,1)
        [trial_corr_mat_this(j,:,:,:),~,shift_mat_this(j,:,:,:)] = trialCorrMat(good_cells,trials(j,:),dat,opt);
        rep_info_this.Stab{j} = squeeze(nanmean(nanmean(trial_corr_mat_this(j,:,1:opt.num_tr_bl,1:opt.num_tr_bl),4),3))';
    end
    corrmat_hi{i} = trial_corr_mat_this;
    shiftmat_hi{i} = shift_mat_this;
    repinfo_hi{i} = rep_info_this;
    
    %% find blocks of trials with gain change + low contrast
    trials = find(dat.trial_gain==opt.gain & dat.trial_contrast==opt.contr_lo);
    trials = reshape(trials,4,numel(trials)/4)';
    trials2 = nan(size(trials,1),size(trials,2)+2*opt.num_tr_bl);
    for j = 1:size(trials,1)
        trials2(j,:) = trials(j,1)-opt.num_tr_bl:trials(j,end)+opt.num_tr_bl;
    end
    trials = trials2(sum(dat.trial_contrast(trials2')~=opt.contr_lo)==0,:); % make sure all correct contrast
    
    % make rep info table for this session
    rep_info_this = table;
    mouse = strsplit(session_name{i},'_');
    rep_info_this.Mouse = repmat(mouse(1),size(trials,1),1);
    rep_info_this.Session = repmat(session_name(i),size(trials,1),1);
    rep_info_this.Rep = (1:size(trials,1))';
    rep_info_this.CellID = repmat({good_cells},size(trials,1),1);
    rep_info_this.Stab = repmat({nan(numel(good_cells),1)},size(trials,1),1);
    
    % compute trial corr mat for each rep
    trial_corr_mat_this = nan(size(trials,1),numel(good_cells),size(trials,2),size(trials,2));
    shift_mat_this = nan(size(trials,1),numel(good_cells),size(trials,2),size(trials,2));
    for j = 1:size(trials,1)
        [trial_corr_mat_this(j,:,:,:),~,shift_mat_this(j,:,:,:)] = trialCorrMat(good_cells,trials(j,:),dat,opt);
        rep_info_this.Stab{j} = squeeze(nanmean(nanmean(trial_corr_mat_this(j,:,1:opt.num_tr_bl,1:opt.num_tr_bl),4),3))';
    end
    corrmat_lo{i} = trial_corr_mat_this;
    shiftmat_lo{i} = shift_mat_this;
    repinfo_lo{i} = rep_info_this;

    %pb.increment();
end
toc;

end