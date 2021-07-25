% Makes the following figures:
% Fig 6B,C
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = paths.intermediate_data;

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% analysis options
opt = load_default_opt;
opt.dark = true;
opt.SpatialBin = 2; % in cm
opt.smoothSigma_dist = 4; % in cm
opt.max_lag = 800; % in cm (for spatial autocorrelation)
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence

% load dist tuning
paths.rez = fullfile(paths.rez,'dark_distance');
binsize = 2;
smoothsig = 4;
dt = load(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',binsize,smoothsig)),'dist_tuning','opt');
dt = dt.dist_tuning;

% which cells to plot
session = 'npI1_0415_dark_1';
cells_to_plot = 61:66;
cells_to_plot = cells_to_plot([5 2 1 6 4 3]); % order them for the fig

%% load data
dat = load(fullfile(paths.data,session));
good_cells = cell_info.CellID(strcmp(cell_info.Session,session));
fr = calcFRVsDist(good_cells,1:max(dat.trial),dat,opt);
fr_zscore = my_zscore(fr);

%% Figure 6B: firing rate snippet
hfig = figure('Position',[200 200 1300 450]);    
hfig.Name = sprintf('Figure 6B: %s c%d-%d firing rate',session,min(cells_to_plot),max(cells_to_plot));
plot_idx = 1+(10000:opt.SpatialBin:16000)/opt.SpatialBin;
xplot = (plot_idx-1)*opt.SpatialBin/100; % (meters)
for cIdx = 1:numel(cells_to_plot)
    subplot(numel(cells_to_plot),1,cIdx);
    y = fr(cells_to_plot(cIdx),plot_idx);
    cIdx_dist = find(strcmp(dt.Session,session) & dt.CellID==good_cells(cells_to_plot(cIdx)));
    if dt.pval(cIdx_dist)<opt.pval_cutoff && dt.prom(cIdx_dist)>opt.min_prom
        plot(xplot,y);
    else
        plot(xplot,y,'k');
    end
    ylim([0 ceil(max(y))]);
    yticks([0 ceil(max(y))])
    box off
    xticks([]);
    set_fig_prefs;
end
xticks(min(xplot):4:max(xplot));
xlabel('Distance run (m)');
ylabel('Firing rate (Hz)');

%% Figure 6C: autocorrelations
hfig = figure('Position',[200 200 400 450]);    
hfig.Name = sprintf('Figure 6C: %s c%d-%d autocorrelation',session,min(cells_to_plot),max(cells_to_plot));
xplot = (-opt.max_lag:opt.SpatialBin:opt.max_lag)/100;
for cIdx = 1:numel(cells_to_plot)
    subplot(numel(cells_to_plot),1,cIdx); hold on;
    y = fr_zscore(cells_to_plot(cIdx),:);
    xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
    cIdx_dist = find(strcmp(dt.Session,session) & dt.CellID==good_cells(cells_to_plot(cIdx)));
    if dt.pval(cIdx_dist)<opt.pval_cutoff && dt.prom(cIdx_dist)>opt.min_prom
        plot(xplot,xc);
    else
        plot(xplot,xc,'k');
    end
    box off
    xticks([]);
    ylim([-0.5 1])
    yticks([0 1])
    plot(xlim,[0 0],'k:')
    plot([0 0],ylim,'k:')
    set_fig_prefs;
end
xticks(min(xplot):4:max(xplot));
xlabel('Lag (m)');
ylabel('Autocorrelation');
