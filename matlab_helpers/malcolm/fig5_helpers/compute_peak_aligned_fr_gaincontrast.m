function [fr_hi,fr_lo] = compute_peak_aligned_fr_gaincontrast(cell_id,reps_trials_hi,reps_trials_lo,dat,opt)
% Computes GC firing rate aligned to peak of baseline FR map for gain
% contrast data
% MGC 3/23/2020

% find blocks of trials with gain change + high contrast
trials_hi = find(dat.trial_gain==opt.gain & dat.trial_contrast==opt.contr_hi);
trials_hi = reshape(trials_hi,4,numel(trials_hi)/4)';
trials2 = nan(size(trials_hi,1),size(trials_hi,2)+2*opt.num_tr_bl);
for j = 1:size(trials_hi,1)
    trials2(j,:) = trials_hi(j,1)-opt.num_tr_bl:trials_hi(j,end)+opt.num_tr_bl;
end
trials_hi = trials2(sum(dat.trial_contrast(trials2')~=opt.contr_hi)==0,:); % make sure all correct contrast

% find blocks of trials with gain change + lo contrast
trials_lo = find(dat.trial_gain==opt.gain & dat.trial_contrast==opt.contr_lo);
trials_lo = reshape(trials_lo,4,numel(trials_lo)/4)';
trials2 = nan(size(trials_lo,1),size(trials_lo,2)+2*opt.num_tr_bl);
for j = 1:size(trials_lo,1)
    trials2(j,:) = trials_lo(j,1)-opt.num_tr_bl:trials_lo(j,end)+opt.num_tr_bl;
end
trials_lo = trials2(sum(dat.trial_contrast(trials2')~=opt.contr_lo)==0,:); % make sure all correct contrast

% high contrast
num_pad = floor(opt.num_bins/2);
fr_hi = nan(size(reps_trials_hi,1),opt.num_bins);
for i = 1:size(reps_trials_hi,1)
    trials_bl = trials_hi(reps_trials_hi(i,1),1:opt.num_tr_bl);
    trial_gc = trials_hi(reps_trials_hi(i,1),opt.num_tr_bl+reps_trials_hi(i,2));
    fr_bl = calcFR(cell_id,trials_bl,dat,opt);
    fr_gc = calcFR(cell_id,trial_gc,dat,opt);
    fr_bl_pad = [nan(1,num_pad) fr_bl nan(1,num_pad)];
    fr_gc_pad = [nan(1,num_pad) fr_gc nan(1,num_pad)];
    [~,maxidx] = nanmax(fr_bl_pad);
    fr_gc_pad_norm = (fr_gc_pad-min(fr_gc_pad))/(max(fr_gc_pad)-min(fr_gc_pad));
    fr_hi(i,:) = fr_gc_pad_norm(maxidx-num_pad:maxidx+num_pad);
end
fr_hi = nanmean(fr_hi,1);

% low contrast
fr_lo = nan(size(reps_trials_lo,1),opt.num_bins);
for i = 1:size(reps_trials_lo,1)
    trials_bl = trials_lo(reps_trials_lo(i,1),1:opt.num_tr_bl);
    trial_gc = trials_lo(reps_trials_lo(i,1),opt.num_tr_bl+reps_trials_lo(i,2));
    fr_bl = calcFR(cell_id,trials_bl,dat,opt);
    fr_gc = calcFR(cell_id,trial_gc,dat,opt);
    fr_bl_pad = [nan(1,num_pad) fr_bl nan(1,num_pad)];
    fr_gc_pad = [nan(1,num_pad) fr_gc nan(1,num_pad)];
    [~,maxidx] = nanmax(fr_bl_pad);
    fr_gc_pad_norm = (fr_gc_pad-min(fr_gc_pad))/(max(fr_gc_pad)-min(fr_gc_pad));
    fr_lo(i,:) = fr_gc_pad_norm(maxidx-num_pad:maxidx+num_pad);
end
fr_lo = nanmean(fr_lo,1);

end

