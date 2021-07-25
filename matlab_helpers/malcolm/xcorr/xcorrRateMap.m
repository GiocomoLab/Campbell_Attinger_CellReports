function rho = xcorrRateMap(cell_id,trials1,trials2,dat,opt)

rateMap1 = calcFR(cell_id,trials1,dat,opt);
rateMap2 = calcFR(cell_id,trials2,dat,opt);

rho = nan(numel(cell_id),opt.max_lag/opt.SpatialBin+1);
for i = 1:numel(cell_id)
    y1 = rateMap1(i,:)-nanmean(rateMap1(i,:));
    y2 = rateMap2(i,:)-nanmean(rateMap2(i,:));
    rho(i,:) = xcorr(y1,y2,opt.max_lag/opt.SpatialBin,'coeff');
end

end