
%% download receptive field data from figshare DOI: 10.25452/figshare.plus.15050289
data_table = readtable('.../receptive_fields/vi_trippy_paris.xlsx')
trippy_path = '.../receptive_fields'
%%
paths = get_paths();

data_path = paths.data;
opt = load_default_opt;
%% grouped by recording on same animal and hemisphere
groups=unique(data_table.group);
marker = {'x','.','o'};
cols = brewermap(3,'Set1');
for iG=1%:numel(groups)
    idx = find(data_table.group == groups(iG));
    animal = data_table.Animal(idx(1));
    probe_nbrs = data_table.Probe_nbr(idx);
    if length(idx)<2
        continue
    end
    
    TC=cell(1,numel(idx));
    RF = cell(1,numel(idx));
    MU = cell(1,numel(idx));
    for ii=1:numel(idx)
        data = load(fullfile(trippy_path,data_table.Trippy_name{idx(ii)},'receptive_fields.mat'));
        
        rf = [];
        rf_noT = [];
        cells_with_rf=[];
        allMU=[];
        for iC=1:numel(data.fields)
            p=data.fields{iC};
            if length(p)>0
                ff= squeeze(data.staMat(:,:,3:7,iC));
                tmp = mean(ff,3);
                tmp = imgaussfilt(tmp);
                Z = zscore(tmp,[],'All');
                rf_noT = cat(3,rf_noT,Z);
                Z(abs(Z)<2.5)=0;
                rf = cat(3,rf,Z);
                cells_with_rf(end+1)=data.good_cells(iC);
                tmpMU=[];
                for iR=1:numel(p)
                    tmpMU=cat(1,tmpMU,p{iR}.mu);
                end
                if size(tmpMU,1)>1
                    tmpMU = mean(tmpMU);
                end
                
                allMU=cat(1,allMU,tmpMU);
                
                
                
            end
            
        end
        spatial_data = load(fullfile(data_path,data_table.Spatial_name{idx(ii)}));
        [corrMat,frMat,~]=trialCorrMat(cells_with_rf,5:20,spatial_data,opt);
        tc = squeeze(mean(frMat,2));
        tcZ = zscore(tc,[],2);
        pdist_tc = pdist(tcZ,'Correlation');
        
        pdist_rf = pdist(reshape(rf_noT,[],size(rf,3))','Correlation');
        v_idx = isfinite(pdist_rf) & isfinite(pdist_tc);
        [a,b]=corrcoef(pdist_rf(v_idx),pdist_tc(v_idx))
        TC{ii}=tc;
        RF{ii}=rf_noT;
        if startsWith(data_table.Hemisphere{idx(ii)},'L')
            allMU(:,1)=176-allMU(:,1);
        end
        MU{ii}=allMU;
    end
    figure('Color','White')
    subplot(1,3,1)
    hold on
    for ii=1:numel(idx)
        scatter(MU{ii}(:,1),MU{ii}(:,2),45,cols(ii,:),'.')
    end
    set(gca,'YDir','reverse')
xlim([0,176])
ylim([0,96])
xlabel('L -> M')
ylabel('D ->U')
title('Individual Receptive Field Centers')
    allMU = cat(1,MU{:});
    allTC = cat(1,TC{:});
    [Y,E]=discretize(allMU(:,1),3);
    
    subplot(1,3,2)
    [~,sid]=sort(allMU(:,1));
    imagesc(zscore(allTC(sid,:),[],2),[-2 2])
    title(groups(iG))
    
    subplot(1,3,3)
    [~,sid]=sort(allMU(:,2));
    imagesc(zscore(allTC(sid,:),[],2),[-2 2])
    title(sprintf('Hemisphere %s',data_table.Hemisphere{idx(1)}))
end
%%
figure('Render','Painters','Color','White','Position',[441   575   735   224])
[~,sid]=sort(allMU(:,1));
    imagesc(1:2:400,1:numel(sid),zscore(allTC(sid,:),[],2),[-2 2])
    colormap(brewermap(20,'BuPu'))
xlabel('Position [cm]')
