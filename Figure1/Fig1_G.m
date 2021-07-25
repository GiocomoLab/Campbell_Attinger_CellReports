% firing rate maps for example cells in small gain manip

% paths = struct;
paths= get_paths;
opt = load_default_opt;

%% examples baseline only
cells_v1 = {
    {'AA4_190804_gain_2.mat',371,5:20},...
    {'AA2_190809_gaincontrast10_2.mat',337,5:20},...
    {'AA1_190728_gaincontrast10_1.mat',191,5:20},...
    {'AA3_190804_gain_1.mat',612,5:20},...
    {'AA1_190728_gaincontrast10_1.mat',215,5:20},...
    {'AA4_190801_gain_1.mat',374,5:20},...
    };

cells_rsc = {
    {'AA2_190807_gain_1.mat',457,5:20},...
    {'AA44_190926_gain_1.mat',305,5:20},...
    {'AA46_190919_gain_1.mat',423,5:20},...
    {'AA47_190927_gaincontrast10_1.mat',406,5:20},...
    {'AA2_190808_gaincontrast10_1.mat',446,5:20},...
    {'AA47_190927_gaincontrast10_1',314,5:20},...
    };

cells_mec = {
    {'npH3_0402_gaincontrast_1',472,5:20},...
    {'npF4_1023_gaincontrast_1',804,5:20},...
    {'npF2_1016_contrasttrack_gainchanges_1',230,53:68},...
    {'npH3_0401_gain_1',849,5:20},...
    {'npI1_0415_gain_1',204,5:20},...
    {'npJ1_0524_gain_1',568,5:20},...
    };



%%
groups = {cells_mec,cells_v1,cells_rsc};

%%
figure('Renderer','Painters','Position',[44         363        1719         442])
for iG = 1:numel(groups)
    cells=groups{iG};
    cntr=0;
    for i = [1 2 7 8 13 14]+(iG-1)*2
        cntr = cntr+1;
        subplot(3,6,i)
        hold on
        trials = cells{cntr}{3};
        dat = load(fullfile(paths.data,cells{cntr}{1}));
        cellidx = find(dat.sp.cids==cells{cntr}{2});
        % compute firing rate and firing rate correlation matrices
        [corrMat,frMat] = trialCorrMat(cells{cntr}{2},trials,dat,opt);
        
        % stability: avg xcorr of first 10 trials
        corrMat = corrMat + diag(nan(size(corrMat,1),1));
        stab_bl = nanmean(corrMat(:));
        
        % make fig
        if isfield(dat.anatomy,'depth')
            depth = dat.anatomy.depth(cellidx);
        else
            
            depth = dat.anatomy.tip_distance(cellidx)-dat.anatomy.FinalDepth+dat.anatomy.z2;
        end
        depth=round(depth);
        if isfield(dat.anatomy,'parent_shifted')
            region = dat.anatomy.parent_shifted{cellidx};
        else
            
            region = dat.anatomy.cluster_parent{cellidx};
        end
        
        celltitle_fig = '';
        spiket = dat.sp.st(dat.sp.clu==cells{cntr}{2});
        [~,~,spike_idx] = histcounts(spiket,dat.post);
  
        imagesc(opt.xbincent,trials,frMat);
        text(200,max(trials)+1.5,sprintf('Stab: %.2f, Peak FR; %d Hz',stab_bl,round(max(frMat,[],'all'))))
        %title(celltitle_fig,'Interpreter','none');
        ylim([min(trials)-.5 max(trials)+.5]);
        yticks([min(trials) max(trials)]);
        ylabel('Trial');
        xticks([1 200 400]);
        xticklabels({'0', 'cm','400'});
        %xlabel('cm');
        xlim([.5 400]);
        
        
       
    end
end
colormap(cbrewer('seq','BuPu',20))
