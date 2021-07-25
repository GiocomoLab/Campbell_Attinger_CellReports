% Makes the following figures:
% Fig S4B,C,D
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% MGC 7/19/2021

%restoredefaultpath;

% add helper functions to path
%addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.sessions = fullfile(paths.intermediate_data,'session_lists');

% load cell_info table
load(fullfile(paths.intermediate_data,fullfile('cell_info','cell_info_Apr2020'))); 

% analysis options
opt = load_default_opt;

session_name = read_session_file('range_of_contrasts_MEC',paths);

% contrast trials
contrast_trials = reshape(11:220,10,21)';
total_trials_per_contrast = 30;

% more params
xbinedges_lick = 0:10:400; % bin edges on track for licking measurements
xbinedges_speed = 0:10:400; % bin edges on track for running speed measurements

% values to compute
lickrate = {};
speed_binned = {};
avg_speed = nan(numel(session_name),numel(opt.contr_all));
num_trials_with_licking = nan(numel(session_name),numel(opt.contr_all)); % proxy for rewards consumed

%% iterate over sessions
for i = 1:numel(session_name)
    % load data
    fprintf('session %d/%d: %s\n',i,numel(session_name),session_name{i});
    load(fullfile(paths.data,session_name{i}));
    
    % running speed vs position (excluding stationary periods)
    speed_this = calcSpeed(posx,opt);
    speed_this_filt = speed_this;  % get rid of stationary periods
    speed_this_filt(speed_this_filt<opt.SpeedCutoff) = nan;
    speed_binned_this = nan(numel(opt.contr_all),numel(xbinedges_speed)-1);
    [~,~,posbin] = histcounts(posx,xbinedges_speed);
    for j = 1:numel(opt.contr_all)
        for k = 1:numel(xbinedges_speed)-1
            speed_binned_this(j,k) = nanmean(speed_this_filt(ismember(trial,contrast_trials) & trial_contrast(trial)==opt.contr_all(j) & posbin==k));
        end
    end
    speed_binned{i} = speed_binned_this;
    
    % average running speed (including stationary periods)
    for j = 1:numel(opt.contr_all)
        total_dist = total_trials_per_contrast * (opt.TrackEnd-opt.TrackStart);
        total_time = sum(ismember(trial,contrast_trials) & trial_contrast(trial)==opt.contr_all(j)) * opt.TimeBin;
        avg_speed(i,j) = total_dist/total_time;
    end
    
    % lick rate
    lickrate_this = nan(numel(opt.contr_all),numel(xbinedges_lick)-1);
    [~,~,lick_idx] = histcounts(lickt,post);
    for j = 1:numel(opt.contr_all)
        keep_licks = ismember(trial(lick_idx),contrast_trials) & trial_contrast(trial(lick_idx))==opt.contr_all(j) & speed_this(lick_idx)>opt.SpeedCutoff; % speed filter removes grooming
        keep_pos = ismember(trial,contrast_trials) & trial_contrast(trial)==opt.contr_all(j) & speed_this>opt.SpeedCutoff;
        lickcounts = histcounts(lickx(keep_licks),xbinedges_lick);
        time_per_bin = histcounts(posx(keep_pos),xbinedges_lick);
        time_per_bin = time_per_bin * opt.TimeBin;
        lickrate_this(j,:) = lickcounts./time_per_bin;       
    end
    lickrate{i} = lickrate_this;
    
    % num trials with licking (proxy for reward consumption)
    % make sure to take the data from the *previous* trial (because licks
    % to consume a trial's reward happen at the beginning of the next
    % trial)
    prev_trial = trial-1;
    prev_trial(prev_trial==0)=1; % fix first trial edge case
    for j = 1:numel(opt.contr_all)
        keep_licks = ismember(prev_trial(lick_idx),contrast_trials) & trial_contrast(prev_trial(lick_idx))==opt.contr_all(j) & speed_this(lick_idx)>opt.SpeedCutoff; % speed filter removes grooming
        num_trials_with_licking(i,j) = numel(unique(trial(lick_idx(keep_licks))));
    end
    
end
percent_rwds = 100 * num_trials_with_licking/total_trials_per_contrast; % percentage of trials on which "reward was consumed" (licking on next trial)

%% average over multiple sessions from the same mouse
mouse = {};
for i = 1:numel(session_name)
    mouse{i} = session_name{i}(3:4);
end
mouse_uniq = unique(mouse);

lickrate_uniq = {};
speed_uniq = {};
percent_rwds_uniq = nan(numel(mouse_uniq),numel(opt.contr_all));
avg_speed_uniq = nan(numel(mouse_uniq),numel(opt.contr_all));
for i = 1:numel(mouse_uniq)
    keep = strcmp(mouse,mouse_uniq{i});
    lickrate_this = lickrate(keep);
    speed_this = speed_binned(keep);
    percent_rwds_this = percent_rwds(keep,:);
    avg_speed_this = avg_speed(keep,:);
    lickrate_tmp = zeros(numel(opt.contr_all),numel(xbinedges_lick)-1);
    speed_tmp = zeros(numel(opt.contr_all),numel(xbinedges_speed)-1);
    percent_rwds_tmp = zeros(1,numel(opt.contr_all));
    avg_speed_tmp = zeros(1,numel(opt.contr_all));
    for j = 1:sum(keep)
        lickrate_tmp = lickrate_tmp + lickrate_this{j};
        speed_tmp = speed_tmp + speed_this{j};
        percent_rwds_tmp = percent_rwds_tmp + percent_rwds_this(j,:);
        avg_speed_tmp = avg_speed_tmp + avg_speed_this(j,:);
    end
    lickrate_uniq{i} = lickrate_tmp/sum(keep);
    speed_uniq{i} = speed_tmp/sum(keep);
    percent_rwds_uniq(i,:) = percent_rwds_tmp/sum(keep);
    avg_speed_uniq(i,:) = avg_speed_tmp/sum(keep);
end
%%
%% lickrate (indiv. mice)
hfigs(1)=figure('Position',[100 100 800 600]);
hfigs(1).Name = 'lickrate indiv mice';
xbincent_lick = xbinedges_lick(1:end-1)+mean(diff(xbinedges_lick))/2;
lickrate_norm_all = nan(numel(opt.contr_all),numel(xbincent_lick),numel(mouse_uniq));
for i = 1:numel(mouse_uniq)
    subplot(4,4,i); hold on;
    lickrate_this = lickrate_uniq{i};
    lickrate_norm = lickrate_this/max(max(lickrate_this));
    lickrate_norm_all(:,:,i) = lickrate_norm;
    for j = 1:numel(opt.contr_all)
        plot(xbincent_lick,lickrate_norm(j,:),'Color',get_color(1,opt.contr_all(j)));
    end
    title(mouse_uniq{i},'Interpreter','none');
    ylabel('norm. lick rate');
    xlabel('cm');
    xlim([200 400]);
    ylim([0 0.5]);
end
%% Figure S4B left: avg norm lick rate vs pos
hfig=figure('Position',[200 200 300 200]); hold on;
hfig.Name = 'Figure S4B left: norm lick rate vs pos';
y = mean(lickrate_norm_all,3);
stderr = std(lickrate_norm_all,[],3)/sqrt(14);
for i = 1:numel(opt.contr_all)
    %errorbar(xbincent_lick,y(i,:),stderr(i,:),'Color',get_color(1,opt.contr_all)(i,:));
    plot(xbincent_lick,y(i,:),'Color',get_color(1,opt.contr_all(i)));
end
xlim([200 400]);
ylim([0 0.3]);
ylabel('norm. lick rate');
xlabel('cm');

% add significance stars
ylim([0 0.35]);
for i = 1:numel(xbincent_speed)
    yy = squeeze(lickrate_norm_all(:,i,:))';
    pval = anova1(yy,[],'off');
    if pval < 0.001
        plot(xbincent_speed(i),max(ylim),'k.');
        plot(xbincent_speed(i),max(ylim)-0.01,'k.');
        plot(xbincent_speed(i),max(ylim)-0.02,'k.');
    elseif pval< 0.01
        plot(xbincent_speed(i),max(ylim),'k.');
        plot(xbincent_speed(i),max(ylim)-0.01,'k.');
    elseif pval < 0.05
        plot(xbincent_speed(i),max(ylim),'k.');
    end
end

%% running speed (indiv. mice)
hfigs(2)=figure('Position',[100 100 800 600]);
hfigs(2).Name = 'running speed indiv mice';
xbincent_speed = xbinedges_speed(1:end-1)+mean(diff(xbinedges_speed))/2;
speed_norm_all = nan(numel(opt.contr_all),numel(xbincent_speed),numel(mouse_uniq));
for i = 1:numel(mouse_uniq)
    subplot(4,4,i); hold on;
    speed_this = speed_uniq{i};
    speed_norm_all(:,:,i) = speed_this/max(max(speed_this));
    for j = 1:numel(opt.contr_all)
        plot(xbincent_speed,speed_this(j,:),'Color',get_color(1,opt.contr_all(j)));
    end
    title(mouse_uniq{i},'Interpreter','none');
    xlim([200 400]);
    ylabel('running speed');
    xlabel('cm');
end

%% Figure S4B right: avg norm speed vs pos
hfig=figure('Position',[200 200 300 200]); hold on;
hfig.Name = 'Figure S4B right: norm speed vs pos';
y = mean(speed_norm_all,3);
stderr = std(speed_norm_all,[],3)/sqrt(14);
for i = 1:numel(opt.contr_all)
    %errorbar(xbincent_speed,y(i,:),stderr(i,:),'Color',get_color(1,opt.contr_all)(i,:));
    plot(xbincent_speed,y(i,:),'Color',get_color(1,opt.contr_all(i)));
end
ylabel('norm. running speed');
xlabel('cm');

% add significance stars
ylim([0.2 0.9]);
for i = 1:numel(xbincent_speed)
    yy = squeeze(speed_norm_all(:,i,:))';
    pval = anova1(yy,[],'off');
    pval = anova1(yy,[],'off');
    if pval < 0.001
        plot(xbincent_speed(i),max(ylim),'k.');
        plot(xbincent_speed(i),max(ylim)-0.02,'k.');
        plot(xbincent_speed(i),max(ylim)-0.04,'k.');
    elseif pval< 0.01
        plot(xbincent_speed(i),max(ylim),'k.');
        plot(xbincent_speed(i),max(ylim)-0.02,'k.');
    elseif pval < 0.05
        plot(xbincent_speed(i),max(ylim),'k.');
    end
end


%% Figure S4C: average running speed by contrast
hfig = figure('Position',[100 100 250 250]); hold on;
hfig.Name = 'Figure S4C: avg speed by contrast';
y = avg_speed_uniq;
means = mean(y);
stderrs = std(y)/sqrt(size(y,1));
for i = 1:numel(mouse_uniq)
    plot(1:numel(opt.contr_all),y(i,:),'k:');
end
errorbar(1:numel(opt.contr_all),means,stderrs,'k-');
xlim([0 numel(opt.contr_all)+1]);
xticks(1:numel(opt.contr_all));
xticklabels(opt.contr_all);
xtickangle(90);
xlabel('contrast');
ylabel('avg. running speed (cm/s)');

%% Figure S4D: percent rewards consumed by contrast
hfig = figure('Position',[100 100 250 250]); hold on;
hfig.Name = 'Figure S4D: pct rewards consumed by contrast';
y = percent_rwds_uniq;
means = mean(y);
stderrs = std(y)/sqrt(size(y,1));
for i = 1:numel(mouse_uniq)
    plot(1:numel(opt.contr_all),y(i,:),'k:');
end
errorbar(1:numel(opt.contr_all),means,stderrs,'k-');
xlim([0 numel(opt.contr_all)+1]);
xticks(1:numel(opt.contr_all));
xticklabels(opt.contr_all);
xtickangle(90);
xlabel('contrast');
ylabel('% rewards consumed');