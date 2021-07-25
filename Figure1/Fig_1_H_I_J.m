paths = get_paths;
paths.figs = fullfile(paths.figs,'fig1');
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
    valid_files{iF}=files(iF).name;
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
        idx = startsWith(reg,'RS');
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
        frMatFlat = zeros(size(frMat,2)*size(frMat,3),size(frMat,1));
        for iT=1:numel(ops.trials)
            idx = (1:200)+(iT-1)*200;
            frMatFlat(idx,:)=squeeze(frMat(:,iT,:))';
        end
        frMatFlatZ = frMatFlat - nanmean(frMatFlat,1);
        xc_temp = zeros(size(corrMat,1),ops.maxLagAutocorr+1);
        for iC=1:size(corrMat,1)
            temp = xcorr(frMatFlatZ(:,iC),ops.maxLagAutocorr,'coeff');
            xc_temp(iC,:)=temp(ops.maxLagAutocorr+1:end);
        end
        
        spMap_this = squeeze(mean(frMat,2))';
        
        for iR=1:numel(regions)
            reg_idx = startsWith(tmp_reg,regions{iR});
            if any(reg_idx)
                MAPS.(regions{iR}) = cat(2,MAPS.(regions{iR}),spMap_this(:,reg_idx));
                DCORR.(regions{iR}) = cat(2,DCORR.(regions{iR}),correlate_distance(:,reg_idx));
                AutoCorr.(regions{iR}) = cat(1,AutoCorr.(regions{iR}),xc_temp(reg_idx,:));
                RATES.(regions{iR}) = cat(2,RATES.(regions{iR}),fr(reg_idx));
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


%% stability distribution, averaged over sites (1H)
cmap=cbrewer('qual','Set2',3,'pchip');
masterfig=figure('Renderer','Painters','Position',[680   810   560   288]);
hold on
pts = 0:0.01:1;
edges = 0:0.05:1;
xi = 0.5*edges(1:end-1)+0.5*edges(2:end);
indfig =figure();
AVG_ALL = [];
for iR = 1:numel(regions)
    figure(indfig)
    subplot(1,3,iR)
    hold on
    uS=unique(FID.(regions{iR}));
    ALL=[];
    for iS=1:numel(uS)
        idx = FID.(regions{iR})==uS(iS);
        if nnz(idx & STAB.(regions{iR})>ops.stab_thresh)>=5
        %[f,xi]=ksdensity(STAB.(regions{iR})(idx),pts);
        a=histcounts(STAB.(regions{iR})(idx),edges);
        a=a/sum(a);
        a(a==0)=nan;
        f=a;
        
        p=nanmean(STAB.(regions{iR})(idx));
        AVG_ALL=cat(1,AVG_ALL,[p,iR]);
        %plot(xi,f,'Color',cmap(iR,:))
        plot(xi,f,'.')
        %pause
        ALL=cat(1,ALL,f);
        end
    end
    figure(masterfig)
    hold on 
    
    boundedline(xi,nanmean(ALL),nanstd(ALL)/sqrt(size(ALL,1)),'cmap',cmap(iR,:),'alpha');
end
xline(0.5,'--','DisplayName','')
xlim([0 1])
xlabel('Stability')
ylabel('Density [cells]')
legend(regions)

figure
anova1(AVG_ALL(:,1),AVG_ALL(:,2))
mec_v=ranksum(AVG_ALL(AVG_ALL(:,2)==1,1),AVG_ALL(AVG_ALL(:,2)==2,1),'tail','left');
mec_r=ranksum(AVG_ALL(AVG_ALL(:,2)==1,1),AVG_ALL(AVG_ALL(:,2)==3,1),'tail','left');
fprintf('MEC: %d, V1: %d, RS: %d \n',nnz(AVG_ALL(:,2)==1),nnz(AVG_ALL(:,2)==2),nnz(AVG_ALL(:,2)==3))
fprintf('MEC vs V1: %.3f \n MEC vs RSC: %.3f \n',mec_v,mec_r)
%% spatial firing rate averaged over sessions (1I)
cmap=cbrewer('qual','Set2',3,'pchip');
figure('Renderer','Painters')
hold on
N=zeros(1,3);
for iR = 1:numel(regions)
    maps = MAPS.(regions{iR});
    fr = RATES.(regions{iR});
    uS=unique(FID.(regions{iR}));
    maps_all=[];
    for iS=1:numel(uS)
        IDX = STAB.(regions{iR})>ops.stab_thresh & FID.(regions{iR})==uS(iS);
    maps_this=maps(:,IDX);
    if nnz(IDX)>5
    %maps = zscore(maps);
    mm=mean(maps_this,2);
    maps_all = cat(2,maps_all,mm);
    end
    end
    mm=mean(maps_all,2);
    err = std(maps_all,[],2)/sqrt(size(maps_all,2));
    boundedline(1:2:400,mm,err,'alpha','cmap',cmap(iR,:))
    N(iR)=size(maps_all,2);
end
for ii=1:4
    xline(80*ii,'--','DisplayName','')
end
legend(regions)
ylim([0 15])


%% autocorr, averaged across sessions (1J)
figure('Position',[ 91   474   847   324],'Renderer','Painters')
subplot(1,2,1)
X=[];
G=[];
x_vec=00:2:ops.maxLagAutocorr*2;
DAT=struct();
peak_idx = ismember(x_vec,[80 160]);
trough_idx = x_vec == 120;
for iR = 1:numel(regions)
    fid = FID.(regions{iR});
    fr = RATES.(regions{iR})';
    stab=STAB.(regions{iR})';
    XC = AutoCorr.(regions{iR})';
    [a,b,c]=unique(fid);
    TMP=nan(ops.maxLagAutocorr+1,numel(a));
    idx_vals = stab>=ops.stab_thresh;
    for iS = 1:numel(a)
        idx_site = c==iS;
        idx = idx_site & idx_vals;
        if nnz(idx)>=5
            TMP(:,iS)=nanmean(XC(:,idx),2);
        end
    end
    mm=nanmean(TMP,2);
    err = nanstd(TMP,[],2)/sqrt(size(TMP,2));
    boundedline(x_vec,mm,err,'alpha','cmap',cmap(iR,:))
    peak = mean(TMP(peak_idx,:));
    trough = TMP(trough_idx,:);
    X=cat(2,X,peak-trough);
    G=cat(2,G,iR*ones(1,numel(peak)));
    DAT.(regions{iR}) = peak-trough;
end
xlim([0,ops.maxLagAutocorr*2])
legend(regions)
box off
subplot(1,2,2)
plotSpread({X(G==1),X(G==2),X(G==3)})
hold on
text(2,-.1,sprintf('%e',ranksum(X(G==1),X(G==2))));
text(3,-.1,sprintf('%e',ranksum(X(G==1),X(G==3))));

