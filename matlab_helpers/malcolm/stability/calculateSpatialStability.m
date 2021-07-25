function rho = calculateSpatialStability(idx,posx,p)
% function to calculate spatial stability by correlating firing rate in
% first half of trial with second half
%
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
% 
% inputs:
%     idx: spike indices
%     posx: position
%     p: params struct (has bin size, track length, etc)
% outputs:
%     rho: spatial stability
    
n = numel(posx);    
fr1 = calculateSmoothedFiringRate(idx(idx<=n/2),posx(1:floor(n/2)),p);
fr2 = calculateSmoothedFiringRate(idx(idx>n/2)-floor(n/2),posx(floor(n/2)+1:n),p);
rho = corr(fr1',fr2');

end