function hfig = make_psth_gaincontrast(cell_id,reps_trials_hi,reps_trials_lo,dat,opt)
%MAKE_PSTH_GAINCONTRAST MGC 3/19/2020

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

% make fig
num_rows = max(size(reps_trials_hi,1),size(reps_trials_lo,1));
xplot_fr = opt.TrackStart:opt.SpatialBin:opt.TrackEnd;
xplot_fr = xplot_fr(1:end-1)+opt.SpatialBin/2;
xplot_xcorr = -opt.max_lag:opt.SpatialBin:opt.max_lag;
hfig = figure('Position',[200 200 900 100 * num_rows]);


% high contrast
ctr = 1;
for i = 1:size(reps_trials_hi,1)
    trials_bl = trials_hi(reps_trials_hi(i,1),1:opt.num_tr_bl);
    trial_gc = trials_hi(reps_trials_hi(i,1),opt.num_tr_bl+reps_trials_hi(i,2));
    fr_bl = calcFR(cell_id,trials_bl,dat,opt);
    fr_gc = calcFR(cell_id,trial_gc,dat,opt);
    xcorr_this = xcorr(fr_gc,fr_bl,opt.max_lag/opt.SpatialBin,'coeff')';
    
    % plot FR
    subplot(num_rows,6,ctr:ctr+1); hold on;
    plot(xplot_fr,fr_bl,'k-');
    plot(xplot_fr,fr_gc,'-','Color',get_color(opt.gain,opt.contr_hi));
    title(sprintf('c=%d, rep%d, trial%d',opt.contr_hi, reps_trials_hi(i,1),reps_trials_hi(i,2)));    
    if i == size(reps_trials_hi,1)
        xlabel('position (cm)');
    else
        xticklabels([]);
    end
    ylabel('FR (Hz)');
    ctr = ctr+2;
    
    % plot xcorr
    subplot(num_rows,6,ctr); hold on;
    plot(xplot_xcorr,xcorr_this,'k.');
    [maxval,maxidx] = max(xcorr_this);
    plot(xplot_xcorr(maxidx),maxval,'ro');
    plot([0 0],ylim(),'k--'); 
    if i == size(reps_trials_hi,1)
        xlabel('lag (cm)');
    end
    ylabel('corr');
    ctr = ctr+4;   
end

% low contrast
ctr = 4;
for i = 1:size(reps_trials_lo,1)
    trials_bl = trials_lo(reps_trials_lo(i,1),1:opt.num_tr_bl);
    trial_gc = trials_lo(reps_trials_lo(i,1),opt.num_tr_bl+reps_trials_lo(i,2));
    fr_bl = calcFR(cell_id,trials_bl,dat,opt);
    fr_gc = calcFR(cell_id,trial_gc,dat,opt);
    xcorr_this = xcorr(fr_gc,fr_bl,opt.max_lag/opt.SpatialBin,'coeff')';
    
    % plot FR
    subplot(num_rows,6,ctr:ctr+1); hold on;
    plot(xplot_fr,fr_bl,'-','Color',get_color(1,opt.contr_lo));
    plot(xplot_fr,fr_gc,'-','Color',get_color(opt.gain,opt.contr_hi));
    title(sprintf('c=%d, rep%d, trial%d',opt.contr_lo,reps_trials_lo(i,1),reps_trials_lo(i,2)));    
    if i == size(reps_trials_lo,1)
        xlabel('position (cm)');
    else
        xticklabels([]);
    end
    ylabel('FR (Hz)');
    ctr = ctr+2;
    
    % plot xcorr
    subplot(num_rows,6,ctr); hold on;
    plot(xplot_xcorr,xcorr_this,'k.');
    [maxval,maxidx] = max(xcorr_this);
    plot(xplot_xcorr(maxidx),maxval,'ro');
    plot([0 0],ylim(),'k--'); 
    if i == size(reps_trials_lo,1)
        xlabel('lag (cm)');
    end
    ylabel('corr');
    ctr = ctr+4;   
end
    

end

