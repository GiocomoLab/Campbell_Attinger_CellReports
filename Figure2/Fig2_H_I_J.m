%% example remapping
opt = load_default_opt;
opt.num_tr_bl = 6;
paths = get_paths;
fig_path = paths.figs;
filenames = {fullfile(paths.data,'npI4_0424_gaincontrast10_smallgainonly_1.mat')};
clusters={[447,389,417, 536,331, 296]};

for iF = 1:numel(filenames)
    trials = 24-6+(0:15);
    dat = load(filenames{iF});
    
    if isfield(dat,'anatomy')
        if isfield(dat.anatomy,'cluster_parent')
            region = dat.anatomy.cluster_parent;
        else
            region = dat.anatomy.parent_shifted;
        end
    else
        continue
    end
    if iscolumn(region)
        region = region';
    end
    

    good_cells = clusters{iF};
    STAB = zeros(size(good_cells));
    if isfield(dat.anatomy,'depth')
        depth  = dat.anatomy.depth(good_idx);
    else
        depth=zeros(size(STAB));
    end
    
    [corrMat,frMatAll]=trialCorrMat(good_cells,trials,dat,opt);
    STAB=nanmean(nanmean(corrMat(:,1:6,1:6),2),3);
    if nnz(STAB>.5)<5
        continue
    end
    % compute firing rate and firing rate correlation matrices
    [~,sid]=sort(STAB,'descend','MissingPlacement','Last');
    [~,sn] = fileparts(filenames{iF});
    % make fig
    n2plot = min(20,size(frMatAll,1));
    hfig = figure('Position',[440   114   526   684],'Renderer','Painters');
    for iC=1:n2plot
        frMat = squeeze(frMatAll(sid(iC),:,:));
        subplot(3,2,iC)
        hold on
        cluID = good_cells(sid(iC));
        rr=region{sid(iC)};
        dd = depth(sid(iC));
        celltitle_save = sprintf('%s_c%d_%s_d=%d_tr%d-%d_stab=%0.2f',sn,cluID,...
            rr,dd,...
            trials(1),trials(end),STAB(sid(iC)));
        celltitle_fig = sprintf('%s c%d\n%s d=%d stab_bl=%0.2f',sn,cluID,...
            rr,dd,STAB(sid(iC)));
        
        
        
        imagesc(opt.xbincent,trials,squeeze(frMat));
        text(50,max(trials)+1.5,sprintf('Stab: %.2f, Peak FR; %d Hz',STAB(sid(iC)),round(max(frMat,[],'all'))))
        ylim([min(trials)-.5 max(trials)+.5]);
        yticks([min(trials) max(trials)]);
        ylabel('Trial');
        xticks([1 200 400]);
        xticklabels({'0', 'cm','400'});
        
        %patches indicating gain value
        for tr = trials
            patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(dat.trial_gain(tr),100),...
                'EdgeColor',get_color(dat.trial_gain(tr),100));
        end
        for tr = 0.5+[trials(opt.num_tr_bl) trials(opt.num_tr_bl+opt.num_tr_gc)]
            plot([0 400],[tr tr],'w-');
        end
        xlim([-15 400]);
        
        
        
        
    end

    

colormap(cbrewer('seq','BuPu',20))
saveas(hfig,fullfile(fig_path,'coherent_shift_ExampleUnits.pdf'))
end
%% calculate trial by trial similarity matrices

loop_tbt_xcorr

%%
savepath = 'tbtxcorr_08';

matfiles = dir(fullfile(savepath,'*.mat'));
ops.max_lag = 20;

regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
SHIFTMat = struct();
Region = struct();
SessionName = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
    SHIFTMat.(regions{iR})=[];
    Region.(regions{iR}) = {};
    SessionName.(regions{iR}) = {};
end


idx_m = triu(true(6),1);
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    [~,sn]=fileparts(matfiles(iF).name);
    
    tmp = data_out.corrMat;
    tmp_bl = nan(size(tmp,1),1);
    for iC=1:numel(tmp_bl)
        tmp_m = squeeze(tmp(iC,1:6,1:6));
        tmp_m = tmp_m+diag(nan(6,1));
        tmp_bl(iC)=nanmean(tmp_m(:));
        
    end
    tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
    stab = [tmp_bl';tmp_bl_gain'];
   

    tmp_reg = data_out.region;
    XC=data_out.corrMat;
    pm_idx = startsWith(tmp_reg,'VISpm');
    agl_idx = startsWith(tmp_reg,'RSPagl');
    tmp_reg(agl_idx)={'aglRSP'};
    tmp_reg(pm_idx)={'VISm'};
    for iR=1:numel(regions)
        reg_idx = startsWith(tmp_reg,regions{iR});
        if ismember(regions{iR},{'MEC','ECT'})
            mult = -1;
        else
            mult = 1;
        end
        if any(reg_idx)
            CorrMat.(regions{iR}) = cat(1,CorrMat.(regions{iR}),XC(reg_idx,:,:));
            FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
            STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
            Region.(regions{iR}) = cat(2,Region.(regions{iR}),tmp_reg(reg_idx));
            SessionName.(regions{iR}) = cat(2,SessionName.(regions{iR}),repmat({sn},1,nnz(reg_idx)));
        end
    end

    
end
%%
ALLC = [];
ALLR = [];
ALLS = [];
ALLSessionName = {};
for iR=1:3
    sites = FID.(regions{iR});
    usites=unique(sites);
    for iS=1:numel(usites)
        sid = sites==usites(iS);
        stab_idx = STAB_BL_Gain.(regions{iR})(1,:)>=.5;
        if nnz(sid & stab_idx)>=5
            idx = sid & stab_idx;
            ALLC=cat(1,ALLC,nanmean(CorrMat.(regions{iR})(idx,:,:)));
            ALLR = cat(1,ALLR,iR);
            ALLSessionName = cat(1,ALLSessionName,unique(SessionName.(regions{iR})(idx)));
        end
    end
end
%%
pca_data=load(fullfile(paths.intermediate_data,'rep_clusters','rep_clusters.mat'));
NAMES = ALLSessionName;
rep_cluster=nan(size(NAMES));
for iN=1:numel(NAMES)
    sn=NAMES{iN};
    sn=strcat(sn(1:end-1),'rep',sn(end));
    NAMES{iN}=sn;
    idx = find(startsWith(pca_data.rep_id,NAMES{iN}),1);
    rep_cluster(iN) = pca_data.rep_cluster(idx);
end
cluster2analyze=1;
%% look at MEC data only for distance
corrMat=CorrMat.MEC;
stab = nanmean(nanmean(corrMat(:,1:6,1:6),2),3)'>=.5;

MEC_clusters= ismember(SessionName.MEC,ALLSessionName(ALLR==1 & rep_cluster~=cluster2analyze));


sid = FID.MEC(MEC_clusters & stab);
names = SessionName.MEC(MEC_clusters & stab);
corrMat  = corrMat(MEC_clusters & stab,:,:);


%% get diagonals and upper halfs of corr matrices
idx = diag(true(1,15),1);
d1=corrMat(:,idx);

upper_half_idx = triu(true(16),1);
upper_half_vecs = corrMat(:,upper_half_idx);


%% distance to all other in site, and all others 2J

idx = triu(true(16),1);
%idx = diag(true(1,15),1);
d1=corrMat(:,idx);
pdist_all = pdist(d1,'corr');
y =squareform(-1*(pdist_all-1));
dd=diag(true(size(y,1),1),0);
y(dd)=nan;
[uS]=unique(sid);

within_distance=nan(1,numel(uS));
across_distance = within_distance;


for iS=1:numel(uS)
    idx = sid==uS(iS);

    this_session = unique(session_names(idx));
    
    
    within_distance(iS)=nanmean(y(idx,idx),'all');
    across_distance(iS) = nanmean(y(idx,~idx),'all');
    %within_distance(iS)=nanmean(within_distance_cell);
    %across_distance(iS)=nanmean(across_distance_cell);
end

figure('Position',[440   317   252   481])


plot([1 2],[within_distance',across_distance'],'k.','MarkerSize',15)
hold on
plot([1 2],[within_distance',across_distance'],'-','Color',[.5 .5 .5])
set(gca,'XTick',[1 2],'XTickLabel',{'within','across'})
xlim([.8 2.2])
box off
title(sprintf('%e',signrank(within_distance,across_distance)));

    %% plot examples 2 H,I
figure('Position',[440   273   250   525],'Renderer','Painters')
trgain = ones(1,16);
trgain(7:10)=.8;
idx = startsWith(names,'npI4_0424_gaincontrast10_smallgainonly_1_1');


 
    this_session = unique(session_names(idx));
    subplot(2,1,1)
    hold on
    f_stab_this=corrMat(idx,:,:);
    pop_mat = squeeze(nanmean(corrMat(idx,:,:)));
    
    imagesc(pop_mat,[0 0.7])
    
    diags_this = d1(idx,:);
    upper_half_this = upper_half_vecs(idx,:);
    nC=nnz(idx);
    cosine_sim  = zeros(nC,1);
    x=pop_mat(upper_half_idx);
    
    for iC=1:nC
        cosine_sim(iC)=dot(x,upper_half_this(iC,:))/(norm(x)*norm(upper_half_this(iC,:)));
        %cosine_sim(iC)=corr(x,upper_half_this(iC,:)');
    end
    title('average similarity')
    [~,sortidx]=sort(cosine_sim);
    axis image
    % patches indicating gain value
    num_tr=16
    for tr = 1:num_tr
        patch([-num_tr/15 0 0 -num_tr/15],[tr tr tr+1 tr+1]-0.5,get_color(trgain(tr),100),...
            'EdgeColor','none');
        patch([tr tr tr+1 tr+1]-0.5,[-num_tr/15 0 0 -num_tr/15],get_color(trgain(tr),100),...
            'EdgeColor','none');
    end
    xlim([-num_tr/15 num_tr]); ylim([-num_tr/15 num_tr]);
    
    % white lines delineated gain change trials
    for tr = 0.5+[opt.num_tr_bl opt.num_tr_bl+opt.num_tr_gc]
        plot([0.5 num_tr+0.5],[tr tr],'w-');
        plot([tr tr],[0.5 num_tr+0.5],'w-');
    end
    
    xlim([-1 num_tr+0.5]);
    ylim([-1 num_tr+0.5]);
    
    subplot(2,1,2)
    imagesc(diags_this(sortidx,:),[0 .7])
    title('first diagonals')
    xlabel(sprintf('within %.2f, across: %.2f, d: %.3f',within_distance(iS),across_distance(iS),within_distance(iS)-across_distance(iS)))
    %pause

box off

    
