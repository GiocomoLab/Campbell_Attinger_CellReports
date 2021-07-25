
paths = get_paths;

opt = load_default_opt;
opt.num_tr_bl = 6;
%% examples with gain change trials
cells_v1 = {

{'AA3_190804_gain_1.mat',612,15:30},...

};


cells_rsc = {

{'AA47_190927_gaincontrast10_1.mat',406,15:30},...
};


cells_mec = {
%{'npF2_1016_contrasttrack_gainchanges_1',230,53:68},...
{'npH3_0401_gain_1',849,15:30},...
{'npF4_1025_gaincontrast_2',845,15:30},...
{'npH3_0402_gaincontrast_1',472,15:30},...

};

    

%%
groups = {cells_mec,cells_v1,cells_rsc};

%%
figure('Renderer','Painters','Position',[44         363        1719         442])
for iG = 1:numel(groups)
    cells=groups{iG};
cntr=1;
for i = [1 2 7 8 13 14]+(iG-1)*2
    
    subplot(3,6,i)
    hold on
    if cntr>size(cells,2)
        continue
    end
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
    ylim([min(trials)-.5 max(trials)+.5]);
    yticks([min(trials) max(trials)]);
    ylabel('Trial');
    xticks([1 200 400]);
    xticklabels({'0', 'cm','400'});
    xlim([.5 400]);

    
    %patches indicating gain value
    for tr = trials
        patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(dat.trial_gain(tr),100),...
            'EdgeColor',get_color(dat.trial_gain(tr),100));
    end
    for tr = 0.5+[trials(opt.num_tr_bl) trials(opt.num_tr_bl+opt.num_tr_gc)]
        plot([0 400],[tr tr],'w-');
    end
    xlim([-15 400]);
    cntr = cntr+1;
end
end
colormap(cbrewer('seq','BuPu',20))
