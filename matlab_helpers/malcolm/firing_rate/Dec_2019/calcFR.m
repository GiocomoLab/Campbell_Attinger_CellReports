function [fr,xbincent] = calcFR(cell_id,trials,dat,opt)
% calculates average firing rate across trials
% MGC 12/16/2019

% SD of gaussian filter for smoothing
smoothSigma = opt.SmoothSigmaFR/opt.SpatialBin;

% extract data for given trials
posx = dat.posx(ismember(dat.trial,trials));
post = dat.post(ismember(dat.trial,trials));

% divide spike counts by occupancy
xbinedges = opt.TrackStart:opt.SpatialBin:opt.TrackEnd;
xbincent = xbinedges(1:end-1)+opt.SpatialBin/2;
time_per_bin = histcounts(posx, xbinedges);
time_per_bin = time_per_bin * opt.TimeBin;

fr = nan(numel(cell_id),numel(xbincent));
for i = 1:numel(cell_id)
    % get spike times for this cell
    spike_t = dat.sp.st(dat.sp.clu==cell_id(i));
    [~,~,spike_idx] = histcounts(spike_t,dat.post);
    
    % only keep spikes in given trials
    spike_t = spike_t(ismember(dat.trial(spike_idx),trials));
    spike_t = spike_t(spike_t>=min(post) & spike_t<=max(post));
    [~,~,spike_idx] = histcounts(spike_t,post);    
    
    % compute firing rate
    fr_this = histcounts(posx(spike_idx), xbinedges);
    fr_this = fr_this./time_per_bin;

    % interpolate missing values
    if sum(isnan(fr_this))>0
        fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
    end

    % smooth firing rate (differently for single trial vs multiple trials)
    numTrials = sum(diff(posx)<-20)+1;
    if numTrials == 1
        % pad beginning and end with same value (instead of zero)
        fr_smooth = [repmat(fr_this(1),1,10) fr_this repmat(fr_this(end),1,10)];
        fr_smooth = gauss_smoothing(fr_smooth,smoothSigma);
        fr_this = fr_smooth(11:end-10);
    else
        % make beginning and end match
        fr_smooth = repmat(fr_this,1,3);
        fr_smooth = gauss_smoothing(fr_smooth,smoothSigma);
        fr_this = fr_smooth(numel(fr_this)+1:numel(fr_this)*2);
    end
    
    fr(i,:) = fr_this;
end

end