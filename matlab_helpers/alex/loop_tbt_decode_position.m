ops = load_default_opt;
ops.trial_range = [-6:9];
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.xbinedges = ops.edges;
ops.xbincent = .5*ops.edges(1:end-1)+.5*ops.edges(2:end);
ops.SpatialBin = ops.BinWidth;
ops.nBins = numel(ops.edges)-1;
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.max_lag = 30;
ops.maxLag = ops.max_lag;
ops.stab_thresh = 0.5;
ops.trials_train = 1:6;
ops.SpeedCutoff = -1;

paths = get_paths();
%%
gain = 0.5;
contrast = 100;
regions = {'MEC','RSP','VISp'};
filenames = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,fullfile(paths.data));
    filenames=cat(2,filenames,tmp1);
    triggers = cat(2,triggers,tmp2);
end

savepath = fullfile('savepath');
if ~isfolder(savepath)
    mkdir(savepath)
end
%

%%
%%
for iF=1:numel(filenames)
    
    try
        [~,sn]=fileparts(filenames{iF});
        
        for iRep=1:numel(triggers{iF})
            data = load(filenames{iF});
            
            if isfield(data.anatomy,'parent_shifted')
                reg = data.anatomy.parent_shifted;
            else
                reg = data.anatomy.cluster_parent;
            end
            if iscolumn(reg)
                reg = reg';
            end
            
            reg=reg(data.sp.cgs==2);
            reg_orig = data.anatomy.cluster_parent((data.sp.cgs==2));
            if iscolumn(reg_orig)
                reg_orig = reg_orig';
            end
            
            
            
            %calculate trial by trial correlation across all bins
            trials=triggers{iF}(iRep)+ops.trial_range;
            ops_here = ops;
            ops_here.trials = trials;
            cellID = data.sp.cids(data.sp.cgs==2);
            reg(startsWith(reg,'VISpm'))={'VISm'};
            reg(startsWith(reg,'RSPa'))={'RA'};
            reg = reg(startsWith(reg,regions));
            cellID = cellID(startsWith(reg,regions));
            
            if numel(cellID)<5
                continue;
            end
            [corrMat,frMat,shiftMat]=trialCorrMat(cellID,trials,data,ops);
            fr = calcFRVsTime(cellID,data,ops);
            speed = calcSpeed(data.posx,ops);
            speed = speed./data.trial_gain(data.trial);
            fr = fr(:,speed>ops.SpeedCutoff);
            trial = data.trial(speed>ops.SpeedCutoff);
            posx = data.posx(speed>ops.SpeedCutoff);
            speed = speed(speed>ops.SpeedCutoff);
            
            % extract firing rate map and position data for trials from this rep
            X = fr(:,ismember(trial,trials));
            mm=max(X,[],2);
            frN=frMat;
            for iT=1:16
                for iP=1:200
                    frN(:,iT,iP)=frMat(:,iT,iP)./mm;
                end
            end
            frMat = frN;
            X=X./mm;
            stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
            if nnz(stab>ops.stab_thresh)<5
                continue
            end
            
            cellID = cellID(stab>ops.stab_thresh);
            region = reg(stab>ops.stab_thresh);
            frMat = frMat(stab>ops.stab_thresh,:,:);
            
            X=X(stab>ops.stab_thresh,:);
            
            Xtilde = X;
            
            
            
            
            corrMat = corrMat(stab>ops.stab_thresh,:,:);
            shiftMat = shiftMat(stab>ops.stab_thresh,:,:);
            score_mat = zeros(2,size(frMat,2),size(frMat,3));
            for iFold = 1:numel(trials)
                take_idx = true(1,numel(trials));
                
                if iFold<=6
                    take_idx(iFold)=false;
                    %calculate tuning curve based on 'training trials'
                end
                take_idx(7:end)=false;
                
                tuning_curve = squeeze(mean(frMat(:,take_idx,:),2));
                vn=vecnorm(tuning_curve,2,1);
                tc_n = tuning_curve./vn;
                Xt=squeeze(frMat(:,iFold,:));
                Xt=Xt./vecnorm(Xt,2,1);
                Xt(isnan(Xt))=0;
                %find closest poin on trajectory
                dot_prod = tc_n' * Xt;
                
                [dist,max_bin] = max(dot_prod);
                tmp_e = ops.xbincent - ops.xbincent(max_bin);
                correction_idx = abs(tmp_e)>ops.TrackEnd/2;
                tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
                score_mat(:,iFold,:)=[tmp_e;dist];
            end
            
            
            trial_this = trial(ismember(trial,trials));
            posx_this = posx(ismember(trial,trials));
            speed_this = speed(ismember(trial,trials));
            posx_shifted = mod(posx_this+200,400);
            [~,~,posbin_orig] = histcounts(posx_this,ops.xbinedges);
            [~,~,posbin] = histcounts(posx_shifted,ops.xbinedges);
            
            % define encoding and decoding trials
            encode_trials = ismember(trial_this,trials(1:ops.num_tr_bl));
            decode_trials = ismember(trial_this,trials(ops.num_tr_bl+1:end));
            Xn = X ./vecnorm(X,2,1);
            dot_prod = tc_n'*Xn;
            [dist,max_bin] = max(dot_prod);
            tmp_e = posx_this' - ops.xbincent(max_bin);
            correction_idx = abs(tmp_e)>ops.TrackEnd/2;
            tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
            

            X = fr(:,ismember(trial,trials));
            X = zscore(X,[],2);
            X=X(stab>ops.stab_thresh,:);
            frMat = zeros(size(X,1),16,200);
            for iT=1:16
                for iPos = 1:200
                    idx = posbin_orig==iPos & trial_this==trials(iT);
                    frMat(:,iT,iPos)=nanmean(X(:,idx),2);
                end
            end
            tc=zeros(size(X,1),200);
            for iPos = 1:200
                idx = posbin_orig==iPos & encode_trials;
                tc(:,iPos)=nanmean(X(:,idx),2);
            end
            
            tuning_curve = squeeze(mean(frMat(:,1:6,:),2))';
            score_mat2 = zeros(size(frMat,2),size(frMat,3));
            for iFold = 1:numel(trials)
                
                Xt=squeeze(frMat(:,iFold,:))';
                [~,iBin]=pdist2(tuning_curve,Xt,'euclidean','Smallest',1);
                tmp_e2 = ops.xbincent - ops.xbincent(iBin);
                correction_idx = abs(tmp_e2)>ops.TrackEnd/2;
                tmp_e2(correction_idx) = tmp_e2(correction_idx)-ops.TrackEnd*sign(tmp_e2(correction_idx));
                score_mat2(iFold,:)=[tmp_e2];
            end
            
            [~,iBin]=pdist2(tuning_curve,X','euclidean','Smallest',1);
            tmp_eGain = posx_this' - ops.xbincent(iBin);
            correction_idx = abs(tmp_eGain)>ops.TrackEnd/2;
            tmp_eGain(correction_idx) = tmp_eGain(correction_idx)-ops.TrackEnd*sign(tmp_eGain(correction_idx));
            
            data_out = matfile(fullfile(savepath,sprintf('%s_%d',sn,iRep)),'Writable',true);
            data_out.corrMat = corrMat;
            data_out.shiftMat = shiftMat;
            data_out.scoreMat = score_mat;
            data_out.scoreMat2 = score_mat2;
            data_out.time_error = tmp_e;
            data_out.time_error2=tmp_eGain;
            %data_out.yhat = yhat;
            %data_out.yhat_error = yhat_error;
            data_out.time_distance = dist;
            data_out.posx = posx_this';
            
            data_out.region = region;
            
            
            data_out.trials = trial_this;
        end
    catch ME
        disp(ME.message)
        disp(sprintf('filenr: %d',iF))
        rethrow(ME)
    end
end