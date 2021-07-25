function [spike_idx_filt,spike_t_filt,posx_filt,post_filt,speed_filt] = ...
    speed_filter(spike_t,posx,post,speed,speed_cutoff)
% function to filter data by running speed (takes out "stationary" periods
% as defined by speed_cutoff)
% MGC 8/23/2019

keep = speed>speed_cutoff;
posx_filt = posx(keep);
post_filt = post(keep);
speed_filt = speed(keep);
[~,~,spike_idx] = histcounts(spike_t,post);
spike_t_filt = spike_t(speed(spike_idx)>speed_cutoff);

% now find indices for speed filtered spikes
spike_idx_filt = nan(size(spike_t_filt));
for i = 1:numel(spike_t_filt)
    [~,spike_idx_filt(i)] = min(abs(spike_t_filt(i)-post_filt));
end

end