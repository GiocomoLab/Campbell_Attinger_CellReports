function speedScoreBinned = speedScoreVR_binned(spike_t,speed,post,posx,p,binsize)
% Computes speed score separately in different spatial bins and averages over bins.
% Bins are weighted by number of spikes in bin.
% Partially deals with correlations between location and running speed.
% MGC 5/30/16
%
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
%
% inputs
%     spike_t: spike times
%     speed: animal's running speed (real coords)
%     post: time stamps
%     posx: positions
%     p: params (time bin, etc)
%     binsize: size of position bins within which to compute speed score
% outputs:
%     speedScoreBinned: speed score

% bin size for speed score
if ~exist('binsize','var')
    binsize = 50;
end

% smoothing parameter for instantaneous firing rate (in time bins)
smoothSigma = 20;
% smoothSigma = 10; % use this to match same smoothing as running speed

% estimate instantaneous firing rate by smoothing spike count histogram
h = hist(spike_t,post);
fr = reshape(h,numel(post),1);
fr = fr./p.TimeBin;
fr = gauss_smoothing(fr,smoothSigma);

% throw out bins with speeds below threshold
select = speed > p.SpeedCutoff; 
speed = reshape(speed,numel(post),1);
speedFilt = speed(select);
frFilt = fr(select);
posxFilt = posx(select);

% bin by position and compute speed score in each position bin
[~,posBin] = histc(posxFilt,p.TrackStart:binsize:p.TrackEnd);
speedScores = nan(max(posBin),1);
for i = 1:max(posBin)
    speedScores(i) = corr(speedFilt(posBin==i),frFilt(posBin==i));
end

% take weighted average across bins, weighted by number of spikes per bin
spikesPerBin = nan(max(posBin),1);
hFilt = h(select);
for i = 1:max(posBin)
    spikesPerBin(i) = sum(hFilt(posBin==i));
end
weights = spikesPerBin/sum(spikesPerBin);
speedScoreBinned = nansum(speedScores.*weights);

end