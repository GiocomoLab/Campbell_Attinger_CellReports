%% run loop_tbt_xcorr with gain=0.5 first
%%
savepath = 'output_of_loop_tbt_xcorr';


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
RepCluster = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    RepCluster.(regions{iR})= [];
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
    
    tmp_sn=strcat(sn(1:end-1),'rep',sn(end));
    idx = find(startsWith(pca_data.rep_id,tmp_sn),1);
    cluster_this = pca_data.rep_cluster(idx);
    if isempty(cluster_this)
        %cluster_this = 0;
        continue
    end
    if ~isfield(data_out,'factors')
        continue
    end
    
    tmp = data_out.corrMat;
    tmp_bl = nan(size(tmp,1),1);
    for iC=1:numel(tmp_bl)
        tmp_m = squeeze(tmp(iC,1:6,1:6));
        tmp_m = tmp_m+diag(nan(6,1));
        tmp_bl(iC)=nanmean(tmp_m(:));
        
    end
    tmp_bl_gain = nanmean(nanmean(tmp(:,7:16,7:16),2),3);
    stab = [tmp_bl';tmp_bl_gain'];
    sm1=data_out.shiftMat; %un altered
    
    
    sm1(isnan(tmp))=nan;
    shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
    shifts = [shift1'];
    %shifts(:,isnan(data_out.factors))=nan;
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
            RepCluster.(regions{iR}) = cat(2,RepCluster.(regions{iR}),cluster_this*ones(1,nnz(reg_idx)));
        end
    end
    
    
end

%%

reg ={'MEC','VISp','RS'}
figure('Position',[441   261   907   538])
for iR=1:3

subplot(1,3,iR)
%scatter(STAB_BL_Gain.(reg{iR})(1,:),STAB_BL_Gain.(reg{iR})(2,:),35,SHIFT.(reg{iR}),'.')
scatter(STAB_BL_Gain.(reg{iR})(1,:),STAB_BL_Gain.(reg{iR})(2,:),35,SHIFT.(reg{iR}),'.')
xlabel('stability BL Pre ')
ylabel('stability gain')
axis image
xlim([-.2 1])
ylim([-.2 1])
set(gca,'CLim',[-3 3])
title(reg{iR})
grid on
colorbar

end

