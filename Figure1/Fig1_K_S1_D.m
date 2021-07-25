%% parameters
ops = load_default_opt;

ops.towerbins = -41:2:41;
ops.towerbincent = ops.towerbins(1:end-1)*.5 + ops.towerbins(2:end)*.5;
ops.xbin = 2;
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
ops.edges = ops.xbinedges;
ops.filter = gausswin(11);
ops.trials=5:20;
ops.num_pcs_decoding = 10;
regions = {'MEC','VISp','RS'};

ops.maxCells = 50; % max number of cells 2 test
ops.nReps = 30; % number of random draws


ops.cells2test = 4:ops.maxCells;

ops.error_edges = 0:4:ops.track_length/2;
ops.error_centers = ops.error_edges(1:end-1)*.5 + ops.error_edges(2:end)*.5;


paths = get_paths;
paths.figs = fullfile(paths.figs,'fig1');
if ~isfolder(paths.figs)
    mkdir(paths.figs)
end

%% load trial by trial spatial firing rate maps, create these maps by running saveSpatialMapsBaseline

files = dir(fullfile(paths.data,'spatialMaps','*.mat'));
FRMAP = [];
Region = {};
for iF=1: numel(files)
    data_out = load(fullfile(files(iF).folder,files(iF).name));
    FRMAP = cat(1,FRMAP,data_out.frMat);
    Region = cat(2,Region,data_out.region);
end
%% allocate matrices for histogram of errors and absolute errors
IDX = struct();

error_histograms = zeros(numel(regions),ops.nReps,numel(ops.error_edges)-1,numel(ops.cells2test));
abs_error = zeros(numel(regions),ops.nReps,numel(ops.cells2test));
position_error = zeros(numel(regions),ops.nReps,numel(ops.xbincent),numel(ops.trials),numel(ops.cells2test));
%% run the loop
cntr = 0;
for nCells2test=ops.cells2test
    cntr = cntr+1;
    
    for iR=1:numel(regions)
        for iRep = 1:ops.nReps
            idx = find(startsWith(Region,regions{iR})); %get index number of all maps for this region
            IDX.(regions{iR})=idx; %save just in case
            
            sub_idx = randsample(idx,nCells2test,false); %sample a random number of these
            X = FRMAP(sub_idx,:,:);
            
            error_pos = [];
            for iFold = 1:numel(ops.trials)
                take_idx = true(1,numel(ops.trials));
                take_idx(iFold)=false;
                %calculate tuning curve based on 'training trials'
                tuning_curve = squeeze(mean(X(:,take_idx,:),2));

                %find closest poin on trajectory
                %dot_prod = tuning_curve' * squeeze(X(:,iFold,:));
                
                %[~,max_bin] = max(dot_prod);
                [~,max_bin]=pdist2(tuning_curve',squeeze(X(:,iFold,:))','euclidean','Smallest',1);
                tmp_e = ops.xbincent - ops.xbincent(max_bin);
                correction_idx = abs(tmp_e)>ops.TrackEnd/2;
                tmp_e(correction_idx) = tmp_e(correction_idx)-ops.TrackEnd*sign(tmp_e(correction_idx));
                error_pos = cat(2,error_pos,tmp_e);
                position_error(iR,iRep,:,iFold,cntr)=tmp_e; %error is already binned, so to average just sum and divide by n trials
            end
            
            
                
            %error_pos(abs(error_pos)>50)=nan;
            N = histcounts(error_pos,ops.error_edges);
            %N = histcounts(error_pos,error_edges);
            error_histograms(iR,iRep,:,cntr)=N;
            abs_error(iR,iRep,cntr) = nanmean(abs(error_pos));
        end
    end
end
%% plot absolute error as function of cells included (1K)
cmap=cbrewer('qual','Set2',3,'pchip');
figure('Color','white')
subplot(10,1,[1:9])
hold on
avgE = squeeze(mean(abs_error,2));
stdE = squeeze(std(abs_error,[],2))./sqrt(size(abs_error,2));
for iR=1:numel(regions)
    boundedline(ops.cells2test,avgE(iR,:),stdE(iR,:),'cmap',cmap(iR,:),'alpha')
end
legend(regions)
xlabel('# of cells')
ylabel('|error|')
box off
set(gca,'TickDir','out');

p_vals = nan(2,47);
for ii=1:47
    p_vals(1,ii)=ranksum(abs_error(1,:,ii),abs_error(2,:,ii));
    p_vals(2,ii)=ranksum(abs_error(2,:,ii),abs_error(3,:,ii));
end
subplot(10,1,10)
imagesc(ops.cells2test,1:2,p_vals>=0.025)
yline(1.5,'w--')
xlim([0 50])
colormap gray

%% plot error histogram (1K inset)
nCells = find(ops.cells2test==20);
normFactor = sum(error_histograms,3);
error_histogramsN=error_histograms;
cntr = 0;
for nCells2test=ops.cells2test
    cntr=cntr+1;
    for iR=1:numel(regions)
        for iRep = 1:ops.nReps
            tmp = squeeze(error_histogramsN(iR,iRep,:,cntr));
            error_histogramsN(iR,iRep,:,cntr)=tmp./sum(tmp);
        end
    end
end

figure('Position',[680   383   560   233],'Renderer','Painters')
subplot(1,2,1)
avg_hist =squeeze(mean((error_histogramsN(:,:,:,nCells)),2));
std_hist = squeeze(std((error_histogramsN(:,:,:,nCells)),[],2))/sqrt(size(error_histogramsN,2));
for iR=1:numel(regions)
    boundedline(ops.error_centers,avg_hist(iR,:),std_hist(iR,:),'cmap',cmap(iR,:),'alpha')
end
legend(regions)
xlabel('|error|')
ylabel('probability')
box off
ylim([0 .05])
set(gca,'TickDir','out');


