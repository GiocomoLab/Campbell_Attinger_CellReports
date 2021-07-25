% this script relies on loop_findshifts first to find shifts
% after, run loop_tbt_xcorr, this will calculate for each site the
% similarity matrices and shift matrices before and after speed correction
savepath = 'tbtxcorr_08'; %output folder of loop_tbt_xcorr

matfiles = dir(fullfile(savepath,'*.mat'));
ops.max_lag = 20;

regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
SHIFTMat = struct();
FACTORS = struct();
Region = struct();
SessionName = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
    SHIFTMat.(regions{iR})=[];
    FACTORS.(regions{iR}) = [];
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
    sm1=data_out.shiftMat; %un altered

    sm1(isnan(tmp))=nan;
    shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
    shifts = [shift1';shift2'];
    shifts(:,isnan(data_out.factors))=nan;
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
            SHIFT.(regions{iR}) = cat(2,SHIFT.(regions{iR}),shifts(:,reg_idx));
            SHIFTMat.(regions{iR}) = cat(1,SHIFTMat.(regions{iR}),sm1(reg_idx,:,:));
            FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
            STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
            FACTORS.(regions{iR}) = cat(2,FACTORS.(regions{iR}),data_out.factors(reg_idx));
            Region.(regions{iR}) = cat(2,Region.(regions{iR}),tmp_reg(reg_idx));
            SessionName.(regions{iR}) = cat(2,SessionName.(regions{iR}),repmat({sn},1,nnz(reg_idx)));
        end
    end
    
    
end

%% get cluster1 sessions
ALLC = [];
ALLSHIFTS = [];
ALLR = [];
ALLS = [];
ALLSessionName = {};
for iR=1:3
    sites = FID.(regions{iR});
    usites=unique(sites);
    for iS=1:numel(usites)
        sid = sites==usites(iS);
        stab_idx = STAB_BL_Gain.(regions{iR})(1,:)>.5;
        if nnz(sid & stab_idx)>=5
            idx = sid & stab_idx;
            ALLC=cat(1,ALLC,nanmean(CorrMat.(regions{iR})(idx,:,:)));
            ALLSHIFTS = cat(2,ALLSHIFTS,nanmean(SHIFT.(regions{iR})(:,idx),2));
            ALLR = cat(1,ALLR,iR);
            ALLS = cat(1,ALLS,nanmean(SHIFTMat.(regions{iR})(idx,:,:)));
            ALLSessionName = cat(1,ALLSessionName,unique(SessionName.(regions{iR})(idx)));
        end
    end
end

%% load clustering data
pca_data=load('rep_clusters.mat');
NAMES = ALLSessionName;
rep_cluster=nan(size(NAMES));
for iN=1:numel(NAMES)
    sn=NAMES{iN};
    sn=strcat(sn(1:end-1),'rep',sn(end));
    NAMES{iN}=sn;
    idx = find(startsWith(pca_data.rep_id,NAMES{iN}),1);
    if ~isempty(idx)
    rep_cluster(iN) = pca_data.rep_cluster(idx);
    end
end
    
%% delta shift
dd={};
X=[];
G=[];
cluster2analyze=1;
for iR=1:3
    tmp = ALLSHIFTS(1,ALLR==iR & rep_cluster == cluster2analyze)-ALLSHIFTS(2,ALLR==iR & rep_cluster == cluster2analyze);
 dd{iR}=tmp;
 X=cat(2,X,tmp);
 G=cat(2,G,iR*ones(size(tmp)));
end
fig = figure();
subplot(1,2,1)
plotSpread(dd)
title('delta shift')
xticklabels(regions)
p=anova1(X,G,'off');



xlabel(sprintf('anova: %.3e',p))
for ii=1:3
    bar(ii,median(X(G==ii)),'EdgeColor','b','FaceColor','none')

    text(ii-.5,1,sprintf('%.3e',signrank(X(G==ii))))
end
text(1,-3,sprintf('%.3e',ranksum(X(G==1),X(G==2))));
text(2,-3,sprintf('%.3e',ranksum(X(G==2),X(G==3))));
box off

p=anova1(X,G);
disp('%%%--- Change in Shift ---%%%')
for ii=1:3
    fprintf('%s diff from 0: %e \n',regions{ii},signrank(X(G==ii)))
end
fprintf('V1 bigger than MEC: %e \n',ranksum(X(G==1),X(G==2)))
fprintf('V1 bigger than RSC: %e \n',ranksum(X(G==2),X(G==3)))
StableClusterNames = ALLSessionName(rep_cluster==cluster2analyze);

%pause
%%
figure(fig)
subplot(1,2,2)
S1=[];
S2=[];
G=[];
for iR=1:3
    tmp1 = ALLSHIFTS(1,ALLR==iR & rep_cluster == cluster2analyze);
    tmp2 = ALLSHIFTS(2,ALLR==iR & rep_cluster == cluster2analyze);
    dd{iR}=tmp2;
    S1 = cat(2,S1,tmp1);
    S2 = cat(2,S2,tmp2);
    G = cat(2,G,iR*ones(size(tmp1)));
end
plotSpread(dd)
p=anova1(S2,G,'off');
xticklabels(regions)
xlabel(sprintf('anova: %.3e',p))
box off
for ii=1:3
    bar(ii,mean(S2(G==ii)),'EdgeColor','b','FaceColor','none')
    text(ii-.5,2,sprintf('%.3e',signrank(S2(G==ii))))
end
title('Residual Shift')
p=anova1(S2,G);
disp('%%%--- Residual Shift ---%%%')
for ii=1:3
    fprintf('%s diff from 0: %e \n',regions{ii},signrank(S2(G==ii)))
end
figure(fig)

text(1,-8,sprintf('MEC: %d,V1:%d,RSC:%d',nnz(G==1),nnz(G==2),nnz(G==3)));
fprintf('V1 smaller than MEC: %e \n',ranksum(S2(G==1),S2(G==2)))
fprintf('RSC smaller than MEC: %e \n',ranksum(S2(G==1),S2(G==3)))
StableClusterNames = ALLSessionName(rep_cluster==cluster2analyze);
saveas(gcf,'/Users/attialex/Dropbox/temporary_images/shift_correction_summary.pdf')