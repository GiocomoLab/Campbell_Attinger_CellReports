function [corrmat,shiftmat,rep_info] = compute_trial_corr_mat_indiv_neurons(session_name,cell_info,paths,opt)
% compute_trial_corr_mat_indiv_neurons.m 
%   MGC 3/30/2020
% remember that this takes max corr over -opt.max_lag:opt.max_lag

corr_mat_all = {};
shift_mat_all = {};
rep_info_all = {};
pb = ParforProgressbar(numel(session_name));
tic
parfor i = 1:numel(session_name)
    % fprintf('analyzing %s session %d/%d: %s\n',opt.brain_region,i,numel(session_name),session_name{i});
    dat = load(fullfile(paths.data,session_name{i}));
    if strcmp(opt.brain_region,'RSP')
        good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{i}) & ...
            contains(cell_info.BrainRegion,opt.brain_region));
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
    corr_mat_all{i} = trial_corr_mat_this;
    shift_mat_all{i} = shift_mat_this;
    rep_info_all{i} = cell_info_this;
    
    pb.increment();
end
toc;

%% concat all cells
corrmat = [];
shiftmat = [];
rep_info = [];
for i = 1:numel(corr_mat_all)
    corrmat = [corrmat; corr_mat_all{i}];
    shiftmat = [shiftmat; shift_mat_all{i}];
    rep_info = [rep_info; rep_info_all{i}];
end

rep_info.rep_id = strcat(rep_info.Session,'_rep',num2str(rep_info.Rep));

end

