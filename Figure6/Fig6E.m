% Makes the following figures:
% Fig 6E
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% 7/19/2021 MGC

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = fullfile(paths.intermediate_data,'dark_distance');

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% load dist tuning
binsize = 2;
smoothsig = 4;
dt = load(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',binsize,smoothsig)),'dist_tuning','opt');
opt = dt.opt;
dt = dt.dist_tuning;

% analysis options
opt.pval_cutoff = 0.01;
opt.min_prom = 0.1; % minimum peak prominence

%% compute autocorr for all dist cells

dist_tuned = dt.pval<opt.pval_cutoff & dt.prom>opt.min_prom;
session = unique(dt.Session);
autocorr = cell(numel(session),1);
depth = cell(numel(session),1);
pb = ParforProgressbar(numel(session));
fprintf('Computing autocorr for %d sessions...\n',numel(session));
tic
parfor sIdx = 1:numel(session)
    dat = load(fullfile(paths.data,session{sIdx}));
    keep = dist_tuned & strcmp(dt.Session,session{sIdx});
    good_cells = dt.CellID(keep);
    fr = calcFRVsDist(good_cells,1:max(dat.trial),dat,opt);
    fr_zscore = my_zscore(fr);
    autocorr_this = nan(numel(good_cells),1+opt.max_lag/opt.SpatialBin);
    depth_this = dt.DepthOnProbe(keep);
    for cIdx = 1:numel(good_cells)
        y = fr_zscore(cIdx,:);
        xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
        xc = xc(opt.max_lag/opt.SpatialBin+1:end);
        autocorr_this(cIdx,:) = xc;
    end   
    autocorr{sIdx} = autocorr_this;
    depth{sIdx} = depth_this;
    pb.increment;
end
toc

%% Figure 6E: this actually plots all 57 sessions, only 4 were chosen for fig
autocorr_x = 0:opt.SpatialBin:opt.max_lag;
for sIdx = 1:numel(session)
    hfig = figure('Position',[300 300 400 500]);
    hfig.Name = session{sIdx};
    [~,sort_idx] = sort(depth{sIdx},'descend');
    y = autocorr{sIdx}(sort_idx,:);
    imagesc(autocorr_x,1:size(y,1),y);
    title(sprintf('%s\n(%d Dist. Cells)',session{sIdx},size(y,1)),'Interpreter','none');
    yticks([]);
    xticks(0:200:800);
    set(gca,'TickDir','out');
    box off
    xlabel('Period (cm)');
    caxis([0 0.4]);
    colorbar;
    colormap summer;
end