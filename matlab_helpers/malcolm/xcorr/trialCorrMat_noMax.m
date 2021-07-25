function [corrMat,frMat] = trialCorrMat_noMax(cell_id,trials,dat,opt)
% same as trialCorrMat but without taking max
% MGC 3/11/2019

max_lag = getOr(opt,'max_lag',0);
num_lags = 1+2*max_lag/opt.SpatialBin;
corrMat = nan(numel(cell_id),numel(trials),numel(trials),num_lags);

frMat = calcTrialFRMat(cell_id,trials,dat,opt); % single trial fr mat
for i = 1:numel(cell_id)
    if numel(cell_id)==1
        fr_this = frMat';
    else
        fr_this = squeeze(frMat(i,:,:))';
    end
    fr_this = fr_this-nanmean(fr_this); % subtract mean on each trial
    xcorr_this = xcorr(fr_this,max_lag/opt.SpatialBin,'coeff')';
    xcorr_this = reshape(xcorr_this,numel(trials),numel(trials),num_lags);    
    corrMat(i,:,:,:) = xcorr_this+repmat(diag(nan(numel(trials),1)),[1 1 num_lags]); % make diagonals nan
end
corrMat = squeeze(corrMat);

end