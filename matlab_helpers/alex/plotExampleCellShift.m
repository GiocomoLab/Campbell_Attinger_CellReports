function plotExampleCellShift(data,cluID,ops,subplot_ax_spikes,subplot_ax_speed)


%%
trials = ops.trials;
nT=numel(trials);
trialMap = nan(1,numel(data.trial_gain));
cntr = 1;
for iT =1:numel(data.trial_gain)
    if ismember(iT,trials)
        trialMap(iT)=cntr;
        cntr=cntr+1;
    end
end
%recreate data.trial map that resets numbers of included trials and sets
%all else to nan
trial_sorted = nan(size(data.trial));
for iT=1:numel(trial_sorted)
    trial_sorted(iT)=trialMap(data.trial(iT));
end

[speed,speed_raw]=calcSpeed(data.posx,ops);
if ~isfield(ops,'speed_filter')
    ops.speed_filter = ops.filter;
end
if ~isempty(ops.speed_filter)
    speed_raw = conv(speed_raw,ops.speed_filter,'same');
end

frMat = calcTrialFRMat(cluID,trials,data,ops); % single trial fr mat

%cellIDX = find(data.sp.cids==cluID);
meanFR = mean(frMat(ops.bl_pre,:));
[tmp_i]=discretize([60,360],ops.edges);
sta=tmp_i(1);
sto=tmp_i(2);
[ma,mi]=max(meanFR(sta:sto));
mi = mi+sta-1;
maxLoc = ops.midpoints(mi);
% extract spike times for this cell
spike_id=data.sp.clu==cluID;
spike_t = data.sp.st(spike_id);
% convert to VR idx
[~,~,spike_idx] = histcounts(spike_t,data.post);
posx=mod(data.posx,max(ops.edges));
col = zeros(numel(spike_idx),3);
for iS = 1:numel(spike_idx)
    cT = data.trial(spike_idx(iS));
    col(iS,:)=get_color(data.trial_gain(cT),data.trial_contrast(cT));
end


ops_temp = ops;
ops_temp.trials = trials;
trial_speed = getSpeedAroundPoint(speed_raw,data.posx,data.trial,ops_temp,maxLoc,ops_temp.speedWindow);
trial_speed_bl_pre = trial_speed(ops.bl_pre);




[sp_blpre,sidx_pre]=sort(trial_speed_bl_pre,'descend');
tmp_rank = 1:numel(ops.bl_pre);
tmp_rank(sidx_pre)=tmp_rank;
tmp_sid = [sidx_pre];

tmp=trial_sorted;
for iT=1:numel(tmp)
    if ~isnan(tmp(iT))
        tmp(iT)=tmp_rank(tmp(iT));
    end
end
xl = [maxLoc-60 maxLoc+40];


if ~isempty(subplot_ax_spikes)
    axes(subplot_ax_spikes)
    scatter(data.posx(spike_idx),tmp(spike_idx),15,col,'.');
    xlim(xl); xline(maxLoc);
    %ylim([0 nT-4])
    set(gca,'YDir','reverse')
    %set(gca,'PlotBoxAspectRatio',[1 .5 1])
    %set(gca,'OuterPosition',[0,0,1,.5])
    xlabel('Position [cm]')
    box off
    axes(subplot_ax_speed)
    scatter(sp_blpre,ops.bl_pre,30,'k','.')
    
    set(gca,'YDir','reverse')
    %set(gca,'OuterPosition',[0,0,1,.5])
    xlabel('Speed [cm/s]')
    set(gca,'PlotBoxAspectRatio',[.5 1 1])
    box off
end

end

