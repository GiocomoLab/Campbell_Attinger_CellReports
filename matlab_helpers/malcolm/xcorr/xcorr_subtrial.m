function [xcorr_all,rep_info] = xcorr_subtrial(session_name, cell_info, opt, paths)
% computes xcorr with sub-trial resolution (splits each trial into chunks)
% pre-baseline: leave one trial out, correlate to rest
% gain change: correlate each gain change trial to all pre-bl trials
% post-baseline: correlate each trial to all pre-bl trials
% stability criteria: average pre-baseline correlation > threshold
% Statistical unit is a rep.
% For a rep to be included, at least a threshold percentage of cells must
% be stable.
% Includes a matched baseline block for comparison.
% 
% 2/3/2019 MGC

%% set up empty variables
max_num_reps = numel(opt.gain) * numel(opt.contrast) * opt.numreps_max * 2; % max num reps in a session
xcorr_all = {};
rep_info = {};

%% iterate over sessions
if datetime(version('-date')) > datetime('Jan 1, 2019')
    ppb = ParforProgressbar(numel(session_name)); % initialize parfor progress bar
end
tic
parfor sesh_idx = 1:numel(session_name)
    %% load data
    fprintf('session %d/%d: %s\n',sesh_idx,numel(session_name),session_name{sesh_idx});
    dat = load(fullfile(paths.data,session_name{sesh_idx}));
    good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name{sesh_idx}));
    
    % empty arrays for xcorr and rep_info for this session
    xcorr_this = nan(max_num_reps,opt.max_num_cells,opt.num_tr_tot*opt.num_pos_chunks,opt.num_lags);
    rep_info_this = table();
    rep_info_this.Session = cell(max_num_reps,1);
    rep_info_this.Gain = nan(max_num_reps,1);
    rep_info_this.Contrast = nan(max_num_reps,1);
    rep_info_this.Rep = nan(max_num_reps,1);
    rep_info_this.GCorBL = cell(max_num_reps,1);
    rep_info_this.NumCells = nan(max_num_reps,1);
    rep_info_this.CellID = nan(max_num_reps,opt.max_num_cells);
    rep_info_this.StabilityBLPre = nan(max_num_reps,opt.max_num_cells);
    rep_info_this.RepID = cell(max_num_reps,1);

    %% iterate over cells
    for cell_idx = 1:numel(good_cells)
        
        % load spike times for this cell
        spike_t = dat.sp.st(dat.sp.clu==good_cells(cell_idx));
        [~,~,spike_idx] = histcounts(spike_t,dat.post);
        
        counter_rep = 1; % counter for reps     
        for gain_idx = 1:numel(opt.gain) % gains
            for contr_idx = 1:numel(opt.contrast) % contrasts                    
                %find trials with this combination of gain + contrast
                gain_change_trials = find(dat.trial_gain==opt.gain(gain_idx) & ...
                    dat.trial_contrast==opt.contrast(contr_idx));
                gain_change_trials = reshape(gain_change_trials,opt.num_tr_gc,...
                    numel(gain_change_trials)/opt.num_tr_gc)';
                numreps_this = min(size(gain_change_trials,1),opt.numreps_max);
                trialblocks = nan(numreps_this,opt.num_tr_tot);
                for block = 1:numreps_this
                    trialblocks(block,:) = gain_change_trials(block,1)-opt.num_tr_bl:...
                        gain_change_trials(block,1)+opt.num_tr_tot-opt.num_tr_bl-1;
                end

                % iterate over reps
                for rep_idx = 1:numreps_this               
                    %% gain change block
                    trialblock = trialblocks(rep_idx,:);
                    baseline_trials = trialblock(1:opt.num_tr_bl);

                    % iterate over trials
                    counter_chunk = 1;
                    for m = 1:numel(trialblock)
                        % single trial firing rate map for this trial
                        tr = trialblock(m);
                        posx_this = dat.posx(dat.trial==tr);
                        post_this = dat.post(dat.trial==tr);
                        spike_t_this = spike_t(spike_t>min(post_this) & spike_t<max(post_this));
                        [~,~,spike_idx_this] = histcounts(spike_t_this,post_this);
                        fr_this = calculateSmoothedFiringRate(spike_idx_this,posx_this,opt);

                        % firing rate map for rest of baseline trials
                        posx_bl = dat.posx(ismember(dat.trial,baseline_trials) & dat.trial~=tr);
                        post_bl = dat.post(ismember(dat.trial,baseline_trials) & dat.trial~=tr);
                        spike_t_bl = spike_t(ismember(dat.trial(spike_idx),baseline_trials) & dat.trial(spike_idx)~=tr);
                        spike_t_bl = spike_t_bl(spike_t_bl>min(post_bl) & spike_t_bl<max(post_bl));
                        [~,~,spike_idx_bl] = histcounts(spike_t_bl,post_bl);
                        fr_bl = calculateSmoothedFiringRate(spike_idx_bl,posx_bl,opt);

                        % within each map chunk, xcorr the two firing rate maps
                        for chunk = 1:opt.num_pos_chunks
                            fr1 = fr_this(opt.pos_chunk_idx(chunk,:));
                            fr2 = fr_bl(opt.pos_chunk_idx(chunk,:));
                            if mean(fr1) > opt.min_fr && mean(fr2) > opt.min_fr % cells must have min fr
                                xcorr_this(counter_rep,cell_idx,counter_chunk,:) = xcorr(fr1-mean(fr1),fr2-mean(fr2),opt.max_lag/opt.SpatialBin,'coeff');
                            end
                            counter_chunk = counter_chunk+1;
                        end
                    end
                    % enter info for this rep
                    rep_info_this.Session{counter_rep} = session_name{sesh_idx};
                    rep_info_this.Gain(counter_rep) = opt.gain(gain_idx);
                    rep_info_this.Contrast(counter_rep) = opt.contrast(contr_idx);
                    rep_info_this.Rep(counter_rep) = rep_idx;
                    rep_info_this.GCorBL{counter_rep} = 'gc';
                    rep_info_this.NumCells(counter_rep) = numel(good_cells);
                    rep_info_this.CellID(counter_rep,cell_idx) = good_cells(cell_idx);
                    rep_info_this.RepID{counter_rep} = ...
                        sprintf('%s_g=%0.1f_c=%d_rep%d_gc',...
                        session_name{sesh_idx},opt.gain(gain_idx),opt.contrast(contr_idx),rep_idx);
                    % enter bl_pre stability
                    rep_info_this.StabilityBLPre(counter_rep,cell_idx) = ...
                        nanmean(squeeze(xcorr_this(counter_rep,cell_idx,1:opt.num_tr_bl*opt.num_pos_chunks,opt.max_lag/opt.SpatialBin+1)));
                    counter_rep = counter_rep+1;

                    %% matched bl block
                    trialblock = trialblocks(rep_idx,:) - opt.num_tr_gc - opt.num_tr_bl;
                    baseline_trials = trialblock(1:opt.num_tr_bl);

                    % double check that these trials really are baseline
                    if all(dat.trial_gain(trialblock)==1 & dat.trial_contrast(trialblock)==opt.contrast(contr_idx))
                        % iterate over trials
                        counter_chunk = 1;
                        for m = 1:numel(trialblock)
                            % single trial firing rate map for this trial
                            tr = trialblock(m);
                            posx_this = dat.posx(dat.trial==tr);
                            post_this = dat.post(dat.trial==tr);
                            spike_t_this = spike_t(spike_t>min(post_this) & spike_t<max(post_this));
                            [~,~,spike_idx_this] = histcounts(spike_t_this,post_this);
                            fr_this = calculateSmoothedFiringRate(spike_idx_this,posx_this,opt);

                            % firing rate map for rest of baseline trials
                            posx_bl = dat.posx(ismember(dat.trial,baseline_trials) & dat.trial~=tr);
                            post_bl = dat.post(ismember(dat.trial,baseline_trials) & dat.trial~=tr);
                            spike_t_bl = spike_t(ismember(dat.trial(spike_idx),baseline_trials) & dat.trial(spike_idx)~=tr);
                            spike_t_bl = spike_t_bl(spike_t_bl>min(post_bl) & spike_t_bl<max(post_bl));
                            [~,~,spike_idx_bl] = histcounts(spike_t_bl,post_bl);
                            fr_bl = calculateSmoothedFiringRate(spike_idx_bl,posx_bl,opt);

                            % within each map chunk, xcorr the two firing rate maps
                            for chunk = 1:opt.num_pos_chunks
                                fr1 = fr_this(opt.pos_chunk_idx(chunk,:));
                                fr2 = fr_bl(opt.pos_chunk_idx(chunk,:));
                                if mean(fr1) > opt.min_fr && mean(fr2) > opt.min_fr % cells must have min fr
                                    xcorr_this(counter_rep,cell_idx,counter_chunk,:) = xcorr(fr1-mean(fr1),fr2-mean(fr2),opt.max_lag/opt.SpatialBin,'coeff');
                                end
                                counter_chunk = counter_chunk+1;
                            end
                        end
                        % enter info for this rep
                        rep_info_this.Session{counter_rep} = session_name{sesh_idx};
                        rep_info_this.Gain(counter_rep) = opt.gain(gain_idx);
                        rep_info_this.Contrast(counter_rep) = opt.contrast(contr_idx);
                        rep_info_this.Rep(counter_rep) = rep_idx;
                        rep_info_this.GCorBL{counter_rep} = 'bl';
                        rep_info_this.NumCells(counter_rep) = numel(good_cells);
                        rep_info_this.CellID(counter_rep,cell_idx) = good_cells(cell_idx);
                        rep_info_this.RepID{counter_rep} = ...
                            sprintf('%s_g=%0.1f_c=%d_rep%d_bl',...
                            session_name{sesh_idx},opt.gain(gain_idx),opt.contrast(contr_idx),rep_idx);
                        % enter bl_pre stability
                        rep_info_this.StabilityBLPre(counter_rep,cell_idx) = ...
                            nanmean(squeeze(xcorr_this(counter_rep,cell_idx,1:opt.num_tr_bl*opt.num_pos_chunks,opt.max_lag/opt.SpatialBin+1)));
                    end
                    counter_rep = counter_rep+1;
                end
            end
        end
    end
    
    % save results to list
    xcorr_all{sesh_idx} = xcorr_this;
    rep_info{sesh_idx} = rep_info_this;

    % increment parfor progressbar
    if datetime(version('-date')) > datetime('Jan 1, 2019')
        ppb.increment();
    end
end
toc
if datetime(version('-date')) > datetime('Jan 1, 2019')
    delete(ppb);
end

% concatenate all sessions
xcorr_all_cat = [];
rep_info_cat = table();
for i = 1:numel(xcorr_all)
    xcorr_all_cat = [xcorr_all_cat; xcorr_all{i}];
    rep_info_cat = [rep_info_cat; rep_info{i}];    
end
xcorr_all = xcorr_all_cat;
rep_info = rep_info_cat;

% filter out empty reps
keep_rep = ~isnan(rep_info.NumCells);
rep_info = rep_info(keep_rep,:);
xcorr_all = xcorr_all(keep_rep,:,:,:);

% keep as few columns of cells as possible
max_num_cells = max(rep_info.NumCells);
xcorr_all = xcorr_all(:,1:max_num_cells,:,:);
rep_info.CellID = rep_info.CellID(:,1:max_num_cells);
rep_info.StabilityBLPre = rep_info.StabilityBLPre(:,1:max_num_cells);