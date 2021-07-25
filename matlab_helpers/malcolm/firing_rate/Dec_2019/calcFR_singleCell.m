function [fr,xbincent] = calcFR_singleCell(posx,spike_idx,opt)
% calculates average firing rate across trials
% MGC 12/16/2019

% SD of gaussian filter for smoothing
smoothSigma = opt.SmoothSigmaFR/opt.SpatialBin;

% divide spike counts by occupancy
xbinedges = opt.TrackStart:opt.SpatialBin:opt.TrackEnd;
xbincent = xbinedges(1:end-1)+opt.SpatialBin/2;
time_per_bin = histcounts(posx, xbinedges);
time_per_bin = time_per_bin * opt.TimeBin;
    
% compute firing rate
fr_this = histcounts(posx(spike_idx), xbinedges);
fr_this = fr_this./time_per_bin;

% interpolate missing values
if sum(isnan(fr_this))>0
fr_this = interp1(find(~isnan(fr_this)),fr_this(~isnan(fr_this)),1:numel(fr_this));
end

% make beginning and end match
fr_smooth = repmat(fr_this,1,3);
fr_smooth = gauss_smoothing(fr_smooth,smoothSigma);
fr = fr_smooth(numel(fr_this)+1:numel(fr_this)*2);

end