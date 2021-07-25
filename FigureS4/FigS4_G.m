
path = 'output_loop_findshifts_gaincontrast';
filenames = dir(fullfile(path,'*gaincontrast*_*.mat'));
regions = {'VISp','RS','MEC'};

ops = load(fullfile(path,'parameters.mat'));
ops = ops.ops;
opt = load_default_opt;

NBlocks = [];
Frac_Good = [];
STABILITY = struct();
GoodCells = struct();
for iF=1:numel(filenames)
    
    filename=filenames(iF).name;
    
    data_out=load(fullfile(path,filename));
    [unique_regions,~,ridx]=unique(data_out.region);
    NBlocks = cat(1,NBlocks,size(data_out.all_stability,1));
    stab_cell = data_out.all_stability;
    stable_blocks = stab_cell>=opt.stab_thresh;
    tmp_factors = data_out.all_factors*-1;
    tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
    
    enough_data_idx = sum(~isnan(tmp_factors(data_out.chunk_contrast==100,:)),1)>=2 & sum(~isnan(tmp_factors(data_out.chunk_contrast==10,:)),1)>=2;%only use data where there are more than 2 stable blocks
    tmp_idx = startsWith(data_out.region,regions);
    stab_cell(~stable_blocks)=nan;
    for iR=1:numel(unique_regions)
        
        r=unique_regions{iR};
        iidx = ridx==iR;
        if startsWith(r,'RSPagl')
            %keyboard
            r = 'RAGL';
        end
        if startsWith(r,'RS')
            r='RSC';
        end
        if startsWith(r,'VISpm')
            r='VISm';
        end
        
        if ~isvarname(r)
            continue
        end
        
        delay_cell_100 = nanmean(tmp_factors(data_out.chunk_contrast==100,iidx))';
        delay_cell_10 = nanmean(tmp_factors(data_out.chunk_contrast==10,iidx))';
        
        %dat(~enough_data_idx(iidx))=nan;
        delay_cell_100(~enough_data_idx(iidx))=nan;
        
        stab_cell_100 = nanmean(data_out.all_stability(data_out.chunk_contrast==100,iidx))';
        stab_cell_10 = nanmean(data_out.all_stability(data_out.chunk_contrast==10,iidx))';
        stab_cell_100(~enough_data_idx(iidx))=nan;
        tmp = [stab_cell_100,stab_cell_10,delay_cell_100,delay_cell_10,ones(size(delay_cell_100))*iF]; %stability, factors,recording_id
        %tmpgc = [sum(~isnan(tmp_factors(:,iidx)),1)', ones(nnz(iidx),1)*size(a.all_stability,1)];
        if ismember(r,fieldnames(STABILITY))
            STABILITY.(r) = cat(1,STABILITY.(r),tmp);
            %GoodCells.(r) = cat(1,GoodCells.(r),tmpgc);
        else
            STABILITY.(r)=tmp;
            %GoodCells.(r)=tmpgc;
        end
    end
    
    
end

%%
figure('Color','White')
regions = {'MEC','RSC','VISp'}
for iR=1:numel(regions)
    [~,~,col] = unique(STABILITY.(regions{iR})(:,5));
    col =ones(size(col));
    subplot(1,3,iR)
    scatter(STABILITY.(regions{iR})(:,3),STABILITY.(regions{iR})(:,4),125,col,'.')
    idx = ~isnan(STABILITY.(regions{iR})(:,3)) & ~isnan(STABILITY.(regions{iR})(:,3));
    aa=unique(STABILITY.(regions{iR})(:,end));
    nrecs= numel(aa);
    valid = nnz(~isnan(STABILITY.(regions{iR})(:,3)) & ~isnan(STABILITY.(regions{iR})(:,3)));
    title(regions{iR})
    xlabel('Delay 100')
    ylabel('Delay 10')
    axis image
    xlim([-.3 .3])
    ylim([-.3 .3])
    grid on
    hold on
    plot([-.3 .3],[-.3 .3],'k--')
    p=signrank(STABILITY.(regions{iR})(:,3),STABILITY.(regions{iR})(:,4));
    text(.1,-.2,sprintf('%.4e \n %d u %d m',p,(valid),nrecs))
    xticks([-.3:0.15:.3])
    yticks([-.3:0.15:.3])
end
saveas(gcf,'/Users/attialex/Dropbox/temporary_images/delay_low_contrast.pdf')
