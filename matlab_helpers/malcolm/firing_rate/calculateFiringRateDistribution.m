function firing_rate_dist = calculateFiringRateDistribution(idx,posx,p)
% Calculates a distribution of firing rates from shuffled spike times
% Shuffles spike times within each trial separately

% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
%
% inputs:
%     idx: spike indices
%     posx: positions
%     p: params (spatial bin size, etc)
% outputs:
%     firing_rate_dist: distribution of shuffled firing rates
%         (bins x shuffles)

% number of shuffles
numShuffles = 1000;

% gaussian filter for smoothing
smoothSigma = p.SmoothSigmaVR/p.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% get trial information
trial = [1 ; cumsum(diff(posx)<-100)+1];
trial_idx = [1 ; find(diff(posx)<-100); numel(trial)];

% compute occupancy
binedges = p.TrackStart:p.SpatialBin:p.TrackEnd;
time_per_bin = histc(posx, binedges);
time_per_bin = time_per_bin * p.TimeBin;

% create distribution of shuffled firing rates
firing_rate_dist = nan(numShuffles,numel(binedges)-1);
for j = 1:numShuffles
    % shuffle spikes within-trial
    idx_shuffle = [];
    for k = 1:numel(trial_idx)-1
        idx_thisTrial = idx(trial(idx)==k);
        trialStartIdx = trial_idx(k)-1;
        trialNumBins = trial_idx(k+1)-trialStartIdx-1;
        idx_shuffle_thisTrial = mod(idx_thisTrial-trialStartIdx+randsample(trialNumBins,1),trialNumBins)+1+trialStartIdx;
        idx_shuffle = [idx_shuffle ; idx_shuffle_thisTrial];
    end
    idx_shuffle = sort(idx_shuffle);
    
    % divide spike counts by occupancy
    firing_rate_shuffle = histcounts(posx(idx_shuffle),binedges)./time_per_bin;

    % smooth shuffled firing rate, padding at ends to avoid taper
    firing_rate_shuffle = conv([repmat(firing_rate_shuffle(1),floor(smoothWindow/2),1); firing_rate_shuffle; ...
        repmat(firing_rate_shuffle(end),floor(smoothWindow/2),1)],gauss_filter,'valid');
    firing_rate_dist(j,:) = firing_rate_shuffle';
end

end
