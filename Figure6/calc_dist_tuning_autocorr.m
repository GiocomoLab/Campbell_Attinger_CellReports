%% Compute distance tuning using autocorrelation
% MGC 4/30/2020

restoredefaultpath;

% add helper functions to path
addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***
paths.rez = fullfile(paths.intermediate_data,'dark_distance');

% load cell_info table
load(fullfile(paths.intermediate_data,'cell_info\cell_info_Apr2020')); 

% get data files
data_files = dir(fullfile(paths.data,'*dark*.mat'));
data_files = {data_files.name}';
data_files = data_files(~contains(data_files,{'with_rewards','AA'}));
mouse_date = cell(size(data_files));
for dIdx = 1:numel(data_files)
    strsplit_this = strsplit(data_files{dIdx},'_');
    mouse_date{dIdx} = sprintf('%s_%s',strsplit_this{1},strsplit_this{2});
end

% analysis options
opt = load_default_opt;
opt.dark = true;
opt.SpatialBin = 2; % in cm
opt.smoothSigma_dist = 4; % in cm
opt.max_lag = 800; % in cm (for spatial autocorrelation)
opt.num_shuf = 300;

%% Iterate over sessions

rez_all = {};
pb = ParforProgressbar(numel(data_files));
fprintf('Computing distance scores for %d sessions...',numel(data_files));
tic
parfor dIdx = 1:numel(data_files)
    session = data_files{dIdx}(1:end-4);
    dat = load(fullfile(paths.data,session));
    good_cells = cell_info.CellID(strcmp(cell_info.Session,session));
    fr = calcFRVsDist(good_cells,1:max(dat.trial),dat,opt);
    fr_zscore = my_zscore(fr);
    
    peak_all = nan(numel(good_cells),1);
    peak_loc_all = nan(numel(good_cells),1);
    peak_prom_all = nan(numel(good_cells),1);
    for cIdx = 1:numel(good_cells)
        y = fr_zscore(cIdx,:);
        xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
        xc = xc(opt.max_lag/opt.SpatialBin+1:end);
        xplot = 0:opt.SpatialBin:opt.max_lag;
        [peaks,locs,~,prominence] = findpeaks(xc);
        if ~isempty(peaks)
            [peak_all(cIdx),idx] = max(peaks);
            peak_loc_all(cIdx) = xplot(locs(idx));
            peak_prom_all(cIdx) = prominence(idx);
        end
    end
    
    peak_shuf = nan(numel(good_cells),opt.num_shuf);
    for sIdx = 1:opt.num_shuf
        fr_shuf = calcFRVsDist_shuf(good_cells,1:max(dat.trial),dat,opt);
        fr_shuf_zscore = my_zscore(fr_shuf);
        for cIdx = 1:numel(good_cells)
            y = fr_shuf_zscore(cIdx,:);
            xc = xcorr(y,y,opt.max_lag/opt.SpatialBin,"coeff");
            xc = xc(opt.max_lag/opt.SpatialBin+1:end);
            peaks = findpeaks(xc);
            if ~isempty(peaks)
                peak_shuf(cIdx,sIdx) = max(peaks);
            end
        end
    end
    
    pval = nan(numel(good_cells),1);
    for cIdx = 1:numel(good_cells)
        pval(cIdx) = sum(peak_shuf(cIdx,:)>=peak_all(cIdx))/opt.num_shuf;
    end
    pval(isnan(peak_all) | sum(isnan(peak_shuf),2)>opt.num_shuf/3) = nan;
    
    rez_this = struct;
    rez_this.peak = peak_all;
    rez_this.period = peak_loc_all;
    rez_this.prom = peak_prom_all;
    rez_this.peak_shuf = peak_shuf;
    rez_this.pval = pval;
    rez_all{dIdx} = rez_this;
    pb.increment;
end
toc

%% Create table of results

dist_tuning = table;
num_row = 300 * numel(rez_all);
dist_tuning.Mouse = cell(num_row,1);
dist_tuning.Date = cell(num_row,1);
dist_tuning.MouseDate = cell(num_row,1);
dist_tuning.Session = cell(num_row,1);
dist_tuning.CellID = nan(num_row,1);
dist_tuning.UniqueID = cell(num_row,1);
dist_tuning.BrainRegion = cell(num_row,1);
dist_tuning.InsertionDepth = nan(num_row,1);
dist_tuning.DepthOnProbe = nan(num_row,1);
dist_tuning.DepthFromSurface = nan(num_row,1);
dist_tuning.peak = nan(num_row,1);
dist_tuning.period = nan(num_row,1);
dist_tuning.prom = nan(num_row,1);
dist_tuning.pval = nan(num_row,1);
dist_tuning.peak_shuf = nan(num_row,opt.num_shuf);
counter = 1;
for dIdx = 1:numel(rez_all)
    session = data_files{dIdx}(1:end-4);
    cIdx = find(strcmp(cell_info.Session,session));
    N = numel(cIdx);
    tIdx = counter:counter+N-1;
    dist_tuning.Mouse(tIdx) = cell_info.Mouse(cIdx);
    dist_tuning.Date(tIdx) = cell_info.Date(cIdx);
    dist_tuning.MouseDate(tIdx) = cell_info.MouseDate(cIdx);
    dist_tuning.Session(tIdx) = cell_info.Session(cIdx);
    dist_tuning.CellID(tIdx) = cell_info.CellID(cIdx);
    dist_tuning.UniqueID(tIdx) = cell_info.UniqueID(cIdx);
    dist_tuning.BrainRegion(tIdx) = cell_info.BrainRegion(cIdx);
    dist_tuning.InsertionDepth(tIdx) = cell_info.InsertionDepth(cIdx);
    dist_tuning.DepthOnProbe(tIdx) = cell_info.DepthOnProbe(cIdx);
    dist_tuning.DepthFromSurface(tIdx) = cell_info.DepthFromSurface(cIdx);
    dist_tuning.peak(tIdx) = rez_all{dIdx}.peak;
    dist_tuning.period(tIdx) = rez_all{dIdx}.period;
    dist_tuning.prom(tIdx) = rez_all{dIdx}.prom;
    dist_tuning.pval(tIdx) = rez_all{dIdx}.pval;
    dist_tuning.peak_shuf(tIdx,:) = rez_all{dIdx}.peak_shuf;
    counter = counter+N;
end
dist_tuning = dist_tuning(1:counter-1,:);

%% save data

save(fullfile(paths.rez,sprintf('dist_tuning_autocorr_bin%d_smooth%d.mat',opt.SpatialBin,opt.smoothSigma_dist)),...
    'dist_tuning','opt');