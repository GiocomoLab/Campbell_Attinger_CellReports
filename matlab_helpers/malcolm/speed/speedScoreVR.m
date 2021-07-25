function [speedScore,params] = speedScoreVR(spike_t,speed,post,p)
% function to calculate speed score in VR using linear fit
% based on Kropff et al 2015
% Malcolm Campbell, 3/30/16
%
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
%
% inputs
%     spike_t: spike times
%     speed: animal's running speed (real coords)
%     post: time stamps
%     p: params (time bin, etc)
% outputs:
%     speedScore: speed score
%     params: [intercept, slope] of linear fit

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

% compute speed score by correlating instantaneous firing rate and speed
speedScore = corr(speedFilt,frFilt);
params = [ones(sum(select),1) speedFilt]\frFilt;

end