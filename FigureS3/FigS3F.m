% Makes the following figures:
% Fig S3F
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

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
opt.br = {'MEC','VISp','RSP'};
opt.br_plot = {'MEC','V1','RSC'}; % brain region names for plotting

% options for decoding
opt.SpeedCutoff = 2;
opt.stab_thresh_decoding = 0.5; % stab thresh for decoding analysis (which cells to keep for decoding)
opt.min_num_stab_cells_decoding = 5;
opt.smoothSigma_time = 0.2;
opt.xbin = 2;
opt.xbinedges = 0:opt.xbin:400;
opt.xbincent = opt.xbinedges(1:end-1)+opt.xbin/2;
opt.max_decoder_error = 50; % cm; remove times when decoder error is greater than this for part of the analysis
opt.distance_metric = 'euclidean';

opt.trialblock_size = 10;

% for regression
opt.xbinedges_regression = 0:20:400;

%% Load sessions to analyze
session_all = cell(numel(opt.br),1);
for i = 1:numel(session_all)
    session_all{i} = read_session_file(sprintf('small_gain_%s',opt.br{i}),paths);
end

%% POSITION DECODING

beta_all = cell(numel(session_all),1);
num_br = numel(opt.br);
for br_idx = 1:num_br
    
    session = session_all{br_idx};
    fprintf('Analyzing %d %s sessions...\n',numel(session),opt.br{br_idx});
    pb = ParforProgressbar(numel(session));
    
    tic
    
    beta_this = nan(numel(session),1);
    parfor session_idx = 1:numel(session)
        
        % load data       
        dat = load(fullfile(paths.data,session{session_idx}));
        
        if strcmp(opt.br{br_idx},'RSP')
            good_cells = cell_info.CellID(strcmp(cell_info.Session,session{session_idx}) & ...
                (strcmp(cell_info.BrainRegion,'RSPd') | strcmp(cell_info.BrainRegion,'RSPv'))); 
        else
            good_cells = cell_info.CellID(strcmp(cell_info.Session,session{session_idx}) & ...
                strcmp(cell_info.BrainRegion,opt.br{br_idx}));
        end

        % calc running speed
        speed = calcSpeed(dat.posx,opt);
        
        % inst fr
        fr = calcFRVsTime(good_cells,dat,opt);
        
        % filter by running speed
        fr = fr(:,speed>opt.SpeedCutoff);
        trial = dat.trial(speed>opt.SpeedCutoff);
        posx = dat.posx(speed>opt.SpeedCutoff);
        speed = speed(speed>opt.SpeedCutoff); 
        
        % find blocks of 10 contiguous high contrast baseline trials
        trialblocks = [];
        bl_trials = dat.trial_gain==1 & dat.trial_contrast==100;
        start_idx = 1;
        keep_going = true;
        while keep_going
            trials_this = start_idx:start_idx+opt.trialblock_size-1;
            if all(bl_trials(trials_this))
                trialblocks = [trialblocks; trials_this];
                start_idx = start_idx + opt.trialblock_size;
            else
                start_idx = start_idx+1;
            end
            if start_idx+opt.trialblock_size-1>numel(bl_trials)
                keep_going = false;
            end
        end
        
        % for regression
        beta_speed = nan(size(trialblocks,1),1);
        pval_speed = nan(size(trialblocks,1),1);
        
        for block_idx = 1:size(trialblocks,1)
            
            trialblock_this = trialblocks(block_idx,:);
            corrmat = trialCorrMat(good_cells,trialblock_this,dat,opt);
            stab = nanmean(nanmean(corrmat,3),2);
            stab_cells = stab>opt.stab_thresh;
            
            if sum(stab_cells)<opt.min_num_stab_cells_decoding
                continue
            else
                % extract firing rate map and position data for trials from this rep
                fr_this = fr(stab_cells,ismember(trial,trialblock_this));
                
                % normalize firing rate
                fr_this = zscore(fr_this,[],2);

                % set nans to 0
                fr_this(isnan(fr_this)) = 0;

                % get trial, position, and speed for this trialblock
                trial_this = trial(ismember(trial,trialblock_this));
                posx_this = posx(ismember(trial,trialblock_this));
                speed_this = speed(ismember(trial,trialblock_this));
                
                [~,~,posbin] = histcounts(posx_this,opt.xbinedges);
                posbin(posbin==0) = 1; 
                
                pred_pos = [];
                true_pos = [];
                speed_decoding = [];
                for tr_idx = 1:opt.trialblock_size
                    % identify decoding and encoding trial
                    decode_trial = ismember(trial_this,trialblock_this(tr_idx));
                    encode_trials = ~decode_trial;

                    % "tuning curve" for encoding trials
                    tc = nan(size(fr_this,1),numel(opt.xbinedges)-1,1); 
                    for i = 1:numel(opt.xbinedges)-1
                        tc(:,i) = mean(fr_this(:,posbin==i & encode_trials),2);
                    end

                    % decode position in decoding trial
                    decode_traj = fr_this(:,decode_trial);
                    euc_dist = pdist2(tc',decode_traj',opt.distance_metric);
                    [~,max_bin] = min(euc_dist);    
                    pred_pos = [pred_pos; opt.xbincent(max_bin)'];
                    true_pos = [true_pos; opt.xbincent(posbin(decode_trial))'];
                    speed_decoding = [speed_decoding; speed_this(decode_trial)];
                end
                
                % circular decoder error
                err = pred_pos-true_pos;
                correction_idx = abs(err)>opt.track_length/2;
                err(correction_idx) = err(correction_idx)-opt.track_length*sign(err(correction_idx));

                % remove times when decoder fails
                err(abs(err)>opt.max_decoder_error) = nan; 
                
                % create regressor matrix
                [~,~,posbin_regress]=histcounts(true_pos,opt.xbinedges_regression);
                X = nan(size(err,1),1+max(posbin_regress));
                X(:,1) = speed_decoding;
                for j = 1:max(posbin_regress)
                    X(:,j+1) = posbin_regress==j;
                end

                % fit linear model, excluding NaNs
                keep = ~isnan(err);
                mdl = fitlm(X(keep,:),err(keep),'Intercept',false);

                beta_speed(block_idx) = mdl.Coefficients.Estimate(1);
                pval_speed(block_idx) = mdl.Coefficients.pValue(1);
            end
        end

        % beta_speed = beta_speed(~isnan(beta_speed));
        % pval_speed = pval_speed(~isnan(pval_speed));
        
        beta_this(session_idx) = nanmean(beta_speed);

        pb.increment();
    end
    toc
    
    beta_all{br_idx} = beta_this;
end


%% *** FLIP SIGN for consistency with sensory delay metric in Fig 4 (otherwise confusing) ***
for i = 1:num_br
    beta_all{i} = -beta_all{i};
end

%% Figure S3F: decoding_error_running_speed_regression_baseline_trials

means = nan(num_br,1);
sems = nan(num_br,1);
for i = 1:num_br
    y = beta_all{i};
    y = y(~isnan(y));
    means(i) = mean(y);
    sems(i) = std(y)/sqrt(numel(y));
end

hfig = figure; hold on;
hfig.Name = 'Figure S3F: decoding_error_running_speed_regression_baseline_trials';
hbar = bar(means);
hbar(1).FaceColor = 'k';
hbar(1).FaceAlpha = 0.2;
errorbar(1:num_br,means,sems,'k.');
for i = 1:num_br
    y = beta_all{i};
    y = y(~isnan(y));
    N = numel(y);
    jit = randn(N,1)*0.05;
    myscatter(i*ones(N,1)+jit,y,'k',0.5);
end
xticks(1:num_br)
xticklabels(opt.br_plot);
ylabel('Regression coefficient (error/speed) (sec)');