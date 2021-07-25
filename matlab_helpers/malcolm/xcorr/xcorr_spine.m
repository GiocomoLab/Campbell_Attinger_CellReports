function corrMat = xcorr_spine(cell_id,trials1,trials2,dat,opt)
% computes cross corr spine between two blocks of trials (trials1 and
% trials2)
% zscore firing rate first
% MGC 4/3/2020

nbins = (opt.TrackEnd-opt.TrackStart)/opt.SpatialBin;
corrMat = nan(numel(cell_id),nbins,nbins);

% compute FR
fr1 = calcFR(cell_id,trials1,dat,opt);
fr2 = calcFR(cell_id,trials2,dat,opt);

% z-score
fr1z = (fr1-mean(fr1,2))./repmat(std(fr1,[],2),1,size(fr1,2));
fr2z = (fr2-mean(fr2,2))./repmat(std(fr2,[],2),1,size(fr2,2));
keyboard;

for i = 1:numel(cell_id)
    corrMat(i,:,:) = fr1z(i,:) * fr2z(i,:)';
end
%corrMat = squeeze(corrMat);

end

