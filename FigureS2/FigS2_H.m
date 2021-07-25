%data_path = '/Volumes/Samsung_T5/attialex/NP_Data/';
paths = get_paths;
data_path = paths.data;

files = dir(fullfile(data_path,'*.mat'));
animal_list = {};
valid_files = {};
for iF=1:numel(files)
    [~,sn]=fileparts(files(iF).name);
    
    parts = strsplit(sn,'_');
    animal = parts{1};
    date = parts{2};
    animal_date = [animal, '_',date];
    valid_files{iF}=fullfile(files(iF).name);
end
%%
ops = load_default_opt;
ops.xbin = 2;
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
lw=ceil(3*ops.smoothSigma_dist);
wx=-lw:lw;
gw=exp(-wx.^2/(2*ops.smoothSigma_dist^2)); gw=gw/sum(gw);
ops.edges = ops.xbinedges;
ops.filter = gw';
ops.trials=5:20;
ops.nTrials2correlate = numel(ops.trials);
ops.maxLagAutocorr = 100;

%%
regions ={'MEC','VISp','RS'};
MAPS = struct();
RATES = struct();
FID = struct();
STAB =struct();
DCORR=struct();
AutoCorr = struct();
N_UNITS = struct();
SPAN = struct();
DURATION = struct();
for iR=1:numel(regions)
    MAPS.(regions{iR})=[];
    RATES.(regions{iR})=[];
    STAB.(regions{iR})= [];
    DCORR.(regions{iR}) = [];
    AutoCorr.(regions{iR}) = [];
    FID.(regions{iR})= [];
    N_UNITS.(regions{iR})=[];
    SPAN.(regions{iR})=[];
    DURATION.(regions{iR})=[];
end
for iF=1:numel(valid_files)
    try
        data = load(fullfile(data_path,valid_files{iF}));
        if ~isfield(data,'anatomy')
            continue
        end
        if isfield(data.anatomy,'parent_shifted')
            reg = data.anatomy.parent_shifted;
        else
            reg = data.anatomy.cluster_parent;
        end
        if iscolumn(reg)
            reg=reg';
        end
        
        idx = startsWith(reg,'VISpm');
        reg(idx)={'VISm'};
        idx = startsWith(reg,'RSPag');
        reg(idx)={'aglRS'};
        idx = startsWith(region,'RS');
        reg(idx)={'RSC'};
        
        
        
        
        good_cells = data.sp.cgs == 2;
        if nnz(good_cells)<2
            continue
        end
        tmp_reg = reg(good_cells);
        good_cells = data.sp.cids(good_cells);
        good_idx = ismember(data.sp.clu,good_cells);
        
        ops.trials = 1:max(data.trial);
        ops.trials = 5:20;
        
        t=range(data.post(ismember(data.trial,ops.trials)));
        start_t = min(data.post(ismember(data.trial,ops.trials)));
        stop_t = max(data.post(ismember(data.trial,ops.trials)));
        t_range = stop_t-start_t;
        fr=zeros(1,numel(good_cells));
        for iC=1:numel(good_cells)
            nSpikes= nnz(data.sp.clu==good_cells(iC) & data.sp.st<stop_t & data.sp.st>start_t);
            fr(iC)=nSpikes/t_range;
        end
        
        [corrMat,frMat,~]=trialCorrMat(good_cells,ops.trials,data,ops);
        bl_idx = data.trial_gain == 1 & data.trial_contrast == 100;
        if numel(bl_idx)<20
            continue
        end
        bl_idx = bl_idx(ops.trials);
        stab_this = nan(1,size(corrMat,1));
        
        correlate_distance = nan(ops.nTrials2correlate,size(corrMat,1));
        
        for iC=1:size(corrMat,1)
            corrmat_this = squeeze(corrMat(iC,:,:));
            
            stab_this(iC)=nanmean(corrmat_this(:));
            for iST=1:ops.nTrials2correlate
                diag_this=diag(corrmat_this,iST);
                correlate_distance(iST,iC)=nanmedian(diag_this);
            end
            
            
        end
        
        
        spMap_this = squeeze(mean(frMat,2))';
        
        for iR=1:numel(regions)
            reg_idx = startsWith(tmp_reg,regions{iR});
            if any(reg_idx)
                MAPS.(regions{iR}) = cat(2,MAPS.(regions{iR}),spMap_this(:,reg_idx));
                DCORR.(regions{iR}) = cat(2,DCORR.(regions{iR}),correlate_distance(:,reg_idx));
                FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
                STAB.(regions{iR}) = cat(2,STAB.(regions{iR}),stab_this(reg_idx));
            end
        end
    catch ME
        rethrow(ME)
        disp(ME.message)
        disp(valid_files{iF})
    end
    
    
end


%% distance correlation S2H
figure('Position',[ 91   474   437   324],'Renderer','Painters')
for iR = 1:numel(regions)
    fid = FID.(regions{iR});
    fr = RATES.(regions{iR})';
    DC = DCORR.(regions{iR});
    [a,b,c]=unique(fid);
    TMP=nan(size(DC,1),numel(a));
    for iS = 1:numel(a)
        TMP(:,iS)=nanmean(DC(:,c==iS & DC(1,:)'>ops.stab_thresh),2);
    end
    TMPN = TMP;%./repmat(TMP(1,:),size(DC,1),1);
    mm=nanmean(TMPN,2);
    err = nanstd(TMP,[],2)/sqrt(size(TMPN,2));
    boundedline(1:size(DC,1),mm,err,'alpha','cmap',cmap(iR,:))
end
xlim([0,size(DC,1)])
legend(regions)
box off

