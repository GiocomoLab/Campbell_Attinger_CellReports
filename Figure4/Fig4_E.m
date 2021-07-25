%% struct_names
paths = get_paths();

%loop_findshifts
path = 'output_of_loop_findshifts';

filenames = dir(fullfile(path,'*.mat'));
ops = load(fullfile(path,'parameters.mat'));
ops = ops.ops;
opt = load_default_opt;
%% files with .8 gain
gain = 0.8;
contrast = 100;
regions = {'VISp','RS','MEC'};
filenames_08 = {};
triggers = {};
for iR = 1:numel(regions)
    
    [tmp1,tmp2] = getFilesCriteria(regions{iR},contrast,gain,paths.data);
    filenames_08=cat(2,filenames_08,tmp1);
    triggers = cat(2,triggers,tmp2);
end

NBlocks = [];
Frac_Good = [];
%%

STABILITY = struct();
GoodCells = struct();
for iF=1:numel(filenames)
    try
        filename=filenames(iF).name;
        if nnz(endsWith(filenames_08,filename))==0
            continue
        end
        a=load(fullfile(path,filename));
        [unique_regions,~,ridx]=unique(a.region);
        NBlocks = cat(1,NBlocks,size(a.all_stability,1));
        stab = a.all_stability;
        stable_blocks = stab>=opt.stab_thresh;
        tmp_factors = a.all_factors;
        tmp_factors(~stable_blocks)=nan; %set shifts to nan where stability <.5
        
        enough_data_idx = sum(~isnan(tmp_factors),1)>=2; %only use data where there are more than 2 stable blocks
        tmp_idx = startsWith(a.region,regions);
        Frac_Good = cat(1,Frac_Good,[sum(~isnan(tmp_factors(:,tmp_idx)),1)', ones(nnz(tmp_idx),1)*size(a.all_stability,1)]);
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
            try
                
                
                dat = nanmean(tmp_factors(:,iidx))';
                dat(~enough_data_idx(iidx))=nan;
                
                stab = nanmean(a.all_stability(:,iidx))';
                tmp = [stab,dat,ones(size(dat))*iF]; %stability, factors,recording_id, firing rate
                tmpgc = [sum(~isnan(tmp_factors(:,iidx)),1)', ones(nnz(iidx),1)*size(a.all_stability,1)];
                if ismember(r,fieldnames(STABILITY))
                    STABILITY.(r) = cat(1,STABILITY.(r),tmp);
                    GoodCells.(r) = cat(1,GoodCells.(r),tmpgc);
                else
                    STABILITY.(r)=tmp;
                    GoodCells.(r)=tmpgc;
                end
            catch ME
                disp(ME.message)
            end
        end
    catch ME
        disp(ME.message)
    end
    
end


%%

reg = {'MEC','VISp','RSC'}
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
AV={};
AV_ALL={};
X=[];
G=[];
figure
for iR = 1:numel(reg)
    data_this = STABILITY.(reg{iR});
    %data_this(data_this<=-.15)=nan;
    [a,b]=unique(data_this(:,3));
    averages  = zeros(size(a));
    for iP = 1:numel(a)
        idx = data_this(:,3)==a(iP);
        averages(iP)=nanmean((data_this(idx,2))*-1);
        X(end+1)=averages(iP);
        G(end+1)=iR;
        if nnz(idx)<nnz(data_this(:,3)==a(iP))*.1
            averages(iP)=nan;
        end
    end
    AV_ALL.(reg{iR}) = averages;
    AV{iR}=averages;
    
end

subplot(1,2,1)
h1 = raincloud_plot(AV{1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
    'box_col_match', 0,'cloud_edge_col', cb(1,:),'line_width',1);
h2 = raincloud_plot(AV{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'cloud_edge_col', cb(4,:),'line_width',1);
h3 = raincloud_plot(AV{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,...
    'cloud_edge_col', cb(5,:),'line_width',1);
legend([h1{1} h2{1} h3{1}], reg);
title('Average shift per recording');
set(gca,'XLim', [-.25 .15], 'YLim', [-20 31]);

box off
set(gcf,'Renderer','Painters')
%saveas(gcf,'F:/temp/figures/shift_raincloud.pdf')
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cb=cb([1 4 5],:);
subplot(1,2,2)
%plotSpread(AV,'distributionColors',cb)
plotSpread(AV)
hold on
%boxplot(X,G,'PlotStyle','compact','Symbol','','Colors',[.5 .5 .5])
z=cellfun(@nanmean,AV);
bar(z,'FaceColor','none','EdgeColor','b')
set(gca,'XTick',[1 2 3],'XTickLabel',reg,'XTickLabelRotation',45)
ylabel('average speed correction lag')
box off
