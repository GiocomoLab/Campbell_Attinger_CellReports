% Figure S6B,C,D
% Identifies and analyzes distance-tuned neurons in tetrode data from the "optic flow track"
% in Campbell et al 2018
% MGC 3/10/2021

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.open_field = fullfile(paths.intermediate_data,'open_field');

opt = struct;
opt.TimeBin = 0.02;
opt.SpatialBin = 2; % in cm
opt.smoothSigma_dist = 4; % in cm
opt.max_lag = 300; % in cm (for spatial autocorrelation)
opt.num_shuf = 300;
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence
opt.grid_score_cutoff = 0.35;
opt.fr_cutoff = 10;

% for plotting example neuron
opt.SmoothSigmaFR = 4;
opt.TrackStart = 0;
opt.TrackEnd = 400;

%% get cells in "optic flow track" condition

cells_all = dir(paths.tetrode_data);
cells_all = {cells_all.name}';
cells_all = cells_all(contains(cells_all,'optic'));

%% get open field data
allCells_OF = readtable(fullfile(paths.open_field,'allCells_OF.csv'));

%% get grid score, border score, and unique cell ID for each neuron
% all neurons were unique

gs = nan(numel(cells_all),1);
bs = nan(numel(cells_all),1);
uid = cell(numel(cells_all),1);
mouse = cell(numel(cells_all),1);
mean_fr = nan(numel(cells_all),1);
for i = 1:numel(cells_all)
    fprintf('getting grid and border score for cell %d/%d: %s\n',i,numel(cells_all),cells_all{i});
    dat = load(fullfile(paths.tetrode_data,cells_all{i}));
    gs(i) = dat.celldata.of_data.grid_score;
    bs(i) = dat.celldata.of_data.border_score;
    uid{i} = dat.celldata.unique_id;
    mouse{i} = dat.celldata.mouse;
    mean_fr(i) = dat.celldata.of_data.mean_rate;
end

%% compute distance tuning metrics (autocorr peak, location, prominence, shuffle peak) for all cells

rng(1); % for consistency

pb = ParforProgressbar(numel(cells_all));
peak_all = nan(numel(cells_all),1);
peak_loc_all = nan(numel(cells_all),1);
peak_prom_all = nan(numel(cells_all),1);
peak_shuf_all = nan(numel(cells_all),opt.num_shuf);
xc_all = nan(numel(cells_all),2*opt.max_lag/opt.SpatialBin+1);
parfor cIdx = 1:numel(cells_all)
    fprintf('Cell %d/%d: %s\n',cIdx,numel(cells_all),cells_all{cIdx});
    dat = load(fullfile(paths.tetrode_data,cells_all{cIdx}));
    
    post = dat.celldata.ws.post;
    posx = dat.celldata.ws.posx;
    total_dist = posx+dat.celldata.track_length*(dat.celldata.ws.trial-1);
    spike_t = dat.celldata.ws.spike_t;

    distbinedges = 0:opt.SpatialBin:max(total_dist);
    distbincent = distbinedges(1:end-1)+opt.SpatialBin/2;

    % time per distance bin
    timeperbin = histcounts(total_dist,distbinedges);
    timeperbin = timeperbin * mean(diff(post));

    % compute distance-binned firing rate
    [~,~,spike_idx] = histcounts(spike_t,post);
    spike_idx = spike_idx(spike_idx>0);
    fr_this = histcounts(total_dist(spike_idx),distbinedges);
    fr_this = fr_this./timeperbin;

    % interpolate missing values
    if sum(isnan(fr_this))>0
        fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
    end

    % smooth firing rate
    fr_this = gauss_smoothing(fr_this,opt.smoothSigma_dist/opt.SpatialBin);

    % compute peak, peak_loc, peak_prom
    y = zscore(fr_this);
    xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
    xc_all(cIdx,:) = xc;
    xc = xc(opt.max_lag/opt.SpatialBin+1:end);
    xplot = 0:opt.SpatialBin:opt.max_lag;
    [peaks,locs,~,prominence] = findpeaks(xc);
    if ~isempty(peaks)
        [peak_all(cIdx),idx] = max(peaks);
        peak_loc_all(cIdx) = xplot(locs(idx));
        peak_prom_all(cIdx) = prominence(idx);
    end
    
    % shuffle
    peak_shuf_this = nan(1,opt.num_shuf);
    for shuf_idx = 1:opt.num_shuf
        spike_t_shuf = mod(spike_t+unifrnd(20,max(post)-20),max(post));
        
        % compute distance-binned firing rate
        [~,~,spike_idx] = histcounts(spike_t_shuf,post);
        spike_idx = spike_idx(spike_idx>0);
        fr_this = histcounts(total_dist(spike_idx),distbinedges);
        fr_this = fr_this./timeperbin;

        % interpolate missing values
        if sum(isnan(fr_this))>0
            fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
        end

        % smooth firing rate
        fr_this = gauss_smoothing(fr_this,opt.smoothSigma_dist/opt.SpatialBin);

        % compute peak, peak_loc, peak_prom
        y = zscore(fr_this);
        xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
        xc = xc(opt.max_lag/opt.SpatialBin+1:end);
        xplot = 0:opt.SpatialBin:opt.max_lag;
        peaks = findpeaks(xc);
        if ~isempty(peaks)
            peak_shuf_this(shuf_idx) = max(peaks);
        end
    end
    peak_shuf_all(cIdx,:) = peak_shuf_this;
    
    pb.increment();
end

%% identify "distance-tuned" cells using same criteria as in paper

dist_pvalue = nan(numel(cells_all),1);
for i = 1:numel(dist_pvalue)
    dist_pvalue(i) = sum(peak_shuf_all(i,:)>=peak_all(i))/opt.num_shuf;
end
dist_cell = dist_pvalue<opt.pval_cutoff & peak_prom_all>opt.min_prom;

% identify grid cells
grid_cell = gs>opt.grid_score_cutoff & mean_fr<opt.fr_cutoff;

% fisher test
contingency_table = crosstab(dist_cell,grid_cell)
[~,fisher_test_pvalue] = fishertest(contingency_table)

%% Figure S6B: example cell
ex_cell = find(grid_cell);
ex_cell = ex_cell(1);

dat = load(fullfile(paths.tetrode_data,cells_all{ex_cell}));
    
post = dat.celldata.ws.post;
posx = dat.celldata.ws.posx;
total_dist = posx+dat.celldata.track_length*(dat.celldata.ws.trial-1);
spike_t = dat.celldata.ws.spike_t;

distbinedges = 0:opt.SpatialBin:max(total_dist);
distbincent = distbinedges(1:end-1)+opt.SpatialBin/2;

% time per distance bin
timeperbin = histcounts(total_dist,distbinedges);
timeperbin = timeperbin * mean(diff(post));

% compute distance-binned firing rate
[~,~,spike_idx] = histcounts(spike_t,post);
spike_idx = spike_idx(spike_idx>0);
fr_this = histcounts(total_dist(spike_idx),distbinedges);
fr_this = fr_this./timeperbin;

% interpolate missing values
if sum(isnan(fr_this))>0
    fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
end

% smooth firing rate
fr_this = gauss_smoothing(fr_this,opt.smoothSigma_dist/opt.SpatialBin);

% compute peak, peak_loc, peak_prom
y = zscore(fr_this);
xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
xc_half = xc(opt.max_lag/opt.SpatialBin+1:end);
xplot = 0:opt.SpatialBin:opt.max_lag;
[peaks,locs,~,prominence] = findpeaks(xc_half);
if ~isempty(peaks)
    [peak_this,idx] = max(peaks);
    peak_loc_this = xplot(locs(idx));
    peak_prom_this = prominence(idx);
end

hfig = figure('Position',[400 400 400 600]);
hfig.Name = sprintf('Figure S6B: Example cell %s',uid{ex_cell}); 

% subplot(3,1,1);
% posx = dat.celldata.ws.posx;
% spike_idx = dat.celldata.ws.spike_idx;
% [fr,xbincent] = calcFR_singleCell(posx,spike_idx,opt);
% plot(xbincent,fr,'r-');
% xticks([0 200 400]);
% xticklabels('');
% set(gca,'FontSize',12);
% box off;
% ylabel('Firing Rate (Hz)');

subplot(3,1,1);
spikepos = posx(spike_idx);
trial = dat.celldata.ws.trial(spike_idx);
myscatter(spikepos,trial,'r',0.2,10);
ylim([0 max(trial)+1]);
box off;
xticks([0 200 400]);
yticks([1 max(trial)]);
xlabel('Track position (cm)');
ylabel('Trial');
set(gca,'FontSize',12);

subplot(3,1,2);
plot(distbincent,fr_this,'r-');
xlim([6200 7200]);
xlabel('Total distance run (cm)');
ylabel('Firing rate (Hz)');
set(gca,'FontSize',12);
box off;

subplot(3,1,3);
x_xcorr = -opt.max_lag:opt.SpatialBin:opt.max_lag;
plot(x_xcorr,xc,'r-');
set(gca,'FontSize',12);
ylabel('Autocorr.');
xlabel('Lag (cm)')
box off;
xlim([-opt.max_lag opt.max_lag]);

%% Figure S6C: plot prominence versus grid score

hfig = figure('Position',[400 400 500 400]); 
hfig.Name = 'Figure S6C: Grid score versus distance peak prominence';
hold on;
myscatter(gs(~grid_cell),peak_prom_all(~grid_cell),'k',0.3,50);
myscatter(gs(grid_cell),peak_prom_all(grid_cell),'r',0.3,50);
xlabel('Grid Score');
ylabel('Distance Peak Prominence');
box off;

keep = ~isnan(peak_prom_all);
[r,p] = corrcoef(gs(keep),peak_prom_all(keep));
title(sprintf('r = %0.2f, p = %0.1e',r(1,2),p(1,2)));

ylim([-0.05 0.6]);
plot([0.35 0.35],[-0.05 0.6],'k--');
set(gca,'FontSize',12);

%% Figure S6D: grid score for distance cells vs non-distance cells

hfig = figure('Position',[200 200 300 500]);
hfig.Name = 'figure S6D: grid score for distance cells vs non-distance cells';
hold on;
myscatter(ones(sum(~dist_cell),1)+randn(sum(~dist_cell),1)*0.02,gs(~dist_cell),'k',0.2);
myscatter(2*ones(sum(dist_cell),1)+randn(sum(dist_cell),1)*0.02,gs(dist_cell),'b',0.2);
plot([0.5 2.5],[opt.grid_score_cutoff opt.grid_score_cutoff],'k--');
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Non-Grid','Grid'})
ylabel('Grid Score')
xticklabels({'Non-Dist','Dist'})
plot([1.1 1.9],[1.3 1.3],'k-')
text(1.5,1.4,'***','HorizontalAlignment','center')
