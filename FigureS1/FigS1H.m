% Makes the following figures:
% Fig S1H
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***

% analysis options
opt = load_default_opt;

session_name = {{'npI4_0425_baseline_1','npI4_0425_dark_with_rewards_1','npI4_0425_baseline_2'},...
    {'npI5_0418_gain_1','npI5_0418_dark_with_rewards_1','npI5_0418_baseline_4'},...
    {'npJ1_0520_baseline_2','npJ1_0520_dark_with_rewards_1','npJ1_0520_baseline_4'}...
    {'npJ1_0521_contrast_1','npJ1_0521_dark_with_rewards_1','npJ1_0521_baseline_1'},...
    {'npJ2_0516_baseline_1','npJ2_0516_dark_with_rewards_1','npJ2_0516_baseline_2'},...
    {'npJ4_0516_baseline_1','npJ4_0516_dark_with_rewards_1','npJ4_0516_baseline_2'}};


% size of blocks in which to compute stability
trialblocksize = 10;

% number of periods
num_periods = numel(session_name{1});

% things to compute
firing_rate_all = {};
stability_all = {};
num_units = nan(numel(session_name),1);

%% iterate over cells
num_trials_all = nan(numel(session_name),num_periods);
for session_num = 1:numel(session_name)
    fprintf('session %d/%d\n',session_num,numel(session_name));
    dat=load(fullfile(paths.data,session_name{session_num}{1}));
    good_cells = dat.sp.cids(dat.sp.cgs==2);
    num_units(session_num) = numel(good_cells);
    firing_rate_this = nan(numel(good_cells),num_periods);
    stability_this = nan(numel(good_cells),num_periods); 
    for j = 1:num_periods
        fprintf('\t%s\n',session_name{session_num}{j});
        if j~=1
            dat = load(fullfile(paths.data,session_name{session_num}{j}));
            good_cells = dat.sp.cids(dat.sp.cgs==2);
        end  
        baseline_trials = find(dat.trial_gain==1 & dat.trial_contrast==100);
        numtrials = floor(numel(baseline_trials)/trialblocksize)*trialblocksize;
        trialblocks = reshape(baseline_trials(1:numtrials),trialblocksize,numtrials/trialblocksize)';
        for k = 1:numel(good_cells)
            spike_t = dat.sp.st(dat.sp.clu==good_cells(k));
            firing_rate_this(k,j) = numel(spike_t)/(max(dat.post)-min(dat.post));
            stability_block = nan(size(trialblocks,1),1);
            for jj = 1:size(trialblocks,1)
                posx_this = dat.posx(ismember(dat.trial,trialblocks(jj,:)));
                post_this = dat.post(ismember(dat.trial,trialblocks(jj,:)));
                spike_t_this = spike_t(spike_t>=min(post_this) & spike_t<max(post_this));
                trial_this = dat.trial(ismember(dat.trial,trialblocks(jj,:)));
                [~,~,idx_this] = histcounts(spike_t_this,post_this);
                stability_block(jj) = calculateSpatialStability_trial(idx_this,posx_this,trial_this-min(trial_this)+1,opt);            
            end
            stability_this(k,j) = nanmean(stability_block);
        end
        num_trials_all(session_num,j) = numtrials;        
    end
    firing_rate_all{session_num} = firing_rate_this;
    stability_all{session_num} = stability_this;
end

%% Figure S1H: mean firing rate and stability for each session
h=figure; 
h.Name = 'FigureS1H: mean firing rate and stability for each session';
subplot(1,2,1); hold on;
means_fr = nan(numel(session_name),num_periods);
for i = 1:numel(session_name)
    means_fr(i,:) = nanmean(firing_rate_all{i});
    plot(1:num_periods,means_fr(i,:),'-o')
end
ylabel('Mean Firing Rate (Hz)');
xlim([0 4]); xticks([1 2 3]); xticklabels({'BL1','DARK','BL2'});
subplot(1,2,2); hold on;
means_stab = nan(numel(session_name),num_periods);
for i = 1:numel(session_name)
    means_stab(i,:) = nanmean(stability_all{i});
    plot(1:num_periods,means_stab(i,:),'-o')
end
ylabel('Mean Spatial Stability');
xlim([0 4]); xticks([1 2 3]); xticklabels({'BL1','DARK','BL2'});
suptitle(sprintf('n=%d sessions (%d units)',numel(session_name),sum(num_units)));
