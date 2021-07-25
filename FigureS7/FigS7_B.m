%% run this first
loop_sliding_window
%% load 05 data
matfiles = dir('output dir of above script');
nChunks = 21;
x_vec=zeros(nChunks,16);
for ii=1:16
    x_vec(:,ii)=linspace(100,300,nChunks)+(ii-1)*400;
    
end
x_vec = reshape(x_vec,1,[]);
CORR_VEC_05=[];
SHIFT_VEC_05=[];
REG=[];
regions={'MEC','VISp','RS'};
ANIMAL={};
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    sn = matfiles(iF).name(1:end-4);
    sn = strcat(sn(1:end-1),'rep',sn(end));
    tmp_idx = (strcmp(pca_data.rep_id,matfiles(iF).name(1:end-4)));
    stab = nanmean(nanmean(data_out.corrMat(:,1:6,1:6),2),3);
    reg_idx = startsWith(data_out.region,regions);
    N=zeros(1,numel(regions));
    for iR=1:numel(regions)
        N(iR)=nnz(startsWith(data_out.region,regions{iR}));
    end
    [~,winner]=max(N);
    
    if isrow(reg_idx)
        reg_idx = reg_idx';
    end
    idx =stab>=.5 & reg_idx;
    
    if nnz(idx)<5
        continue
    end
    nC=numel(stab);
    nT = size(data_out.corrMat,2);
    nChunks = size(data_out.PEAKS,3);
    
    corr_vec =zeros(nC,nT*nChunks);
    shift_vec=corr_vec;
    PEAKS = data_out.PEAKS;
    SHIFTS = data_out.SHIFTS;
    %valid_idx = ismember(SHIFTS,[-30 30]);
    %SHIFTS(valid_idx)=nan;
    %PEAKS(valid_idx)=nan;
    
    for iC=1:nC
        tmp = reshape(squeeze(PEAKS(iC,:,:))',1,[]);
        corr_vec(iC,:)=tmp;
        tmp = reshape(squeeze(SHIFTS(iC,:,:))',1,[]);
        shift_vec(iC,:)=tmp;
        
    end
    
    CORR_VEC_05 = cat(1,CORR_VEC_05,nanmean(corr_vec(idx,:)));
    SHIFT_VEC_05 = cat(1,SHIFT_VEC_05,nanmean(shift_vec(idx,:)));
    REG = cat(1,REG,winner);
    parts=strsplit(sn,'_');
    ANIMAL = cat(1,ANIMAL,parts{1});
end


%% only take non overlapping chunks
nChunks = 21;
sel_idx=false(nChunks,16);
sel_idx(1,:)=true;
sel_idx(21,:)=true;
sel_idx = sel_idx(:);



CT=cbrewer('qual', 'Set2', 3);
figure('Color','white','Renderer','Painters')
for ii=1:3
    reg_idx = REG==ii;
    subplot(2,1,1)
    hold on
    tmp_corr = mean(CORR_VEC_05(reg_idx,:));
    tmp_shift = mean(SHIFT_VEC_05(reg_idx,:));
    boundedline(x_vec(sel_idx),tmp_corr(sel_idx),std(CORR_VEC_05(reg_idx,(sel_idx)))/sqrt(nnz(reg_idx)),'cmap',CT(ii,:),'alpha')
    
    %plot(x_vec,tmp_corr,'Color',CT(iS+2,:))
    subplot(2,1,2)
    hold on
    boundedline(x_vec(sel_idx),tmp_shift(sel_idx),std(SHIFT_VEC_05(reg_idx,(sel_idx)))/sqrt(nnz(reg_idx)),'cmap',CT(ii,:),'alpha')
    tmp = CORR_VEC_05(reg_idx,1);
    
    uA = unique(ANIMAL(reg_idx & ~isnan(CORR_VEC_05(:,1))));
    fprintf('%s: %d blocks, %d mice \n',reg{ii},nnz(~isnan(tmp)),numel(uA));
end

for iS=1:2
    subplot(2,1,iS)
    for ii=1:4
        xline((5+ii)*400+1);
    end
end
subplot(2,1,2)
title('shift')
% boundedline(x_vec,mean(SHIFT_VEC_05),std(SHIFT_VEC_05)/sqrt(size(SHIFT_VEC_05,1)),'cmap',[.4 .4 .4])
subplot(2,1,1)
title('similarity')
ax1=gca;
ax1.XAxis.TickLength = [0 0];
ax1.XAxis.Limits=[0 16*400];
legend(regions)

