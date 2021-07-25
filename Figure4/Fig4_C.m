%% this script requires that shifts are estimated first with 
paths = get_paths();
data_path = paths.data;
if ~isfolder(save_path)
    mkdir(save_path);
end
ops = load(fullfile(shift_path,'parameters.mat'));
ops = ops.ops;
ops.BinWidth = 2;
ops.SpatialBin=2;
ops.edges = 0:ops.BinWidth:400;
ops.nBins = numel(ops.edges)-1;

ops.speedWindow = [-10 -1]; % in cm
ops.smoothSigma=ops.smoothSigma_dist;
smoothSigma = ops.smoothSigma/ops.BinWidth;
ops.filter = gausswin(floor(smoothSigma*5/2)*2+1);
ops.filter = ops.filter/sum(ops.filter);
ops.midpoints = ops.edges(1:end-1)*.5+ops.edges(2:end)*.5;

%%
cells = {
{'npF4_1023_gaincontrast_1.mat',804,6:20};...%MEC
{'AA1_190726_gain_1.mat',380,102:116};... %V1
{'AA44_190927_gaincontrast10_1.mat',131,6:20};...%RSC
{'npI4_0424_gaincontrast10_smallgainonly_1.mat',541,9:23};... %MEC
{'AA46_190921_gaincontrast10_combined.mat',521,6:20};...$V1
{'AA47_190927_gaincontrast10_1.mat',301,6:20};...%RSC
{'AA4_190805_gaincontrast10_2.mat',321,6:20};...%V1 example shift
};

%%
masterfig = figure('Position',[240         791        1518         307]);

for iC=1:numel(cells)
    data = load(fullfile(data_path,cells{iC}{1}));
    
    cluID = cells{iC}{2};
    ops.trials = cells{iC}{3};
    ops.bl_pre = 1:15;
        ops.gain_trials = 16:19;
        ops.bl_post = [];
     if iC<=6
     figure(masterfig)
     subplot_ax_spikes = subplot(2,6,iC*2-1);
     subplot_ax_speed = subplot(2,6,iC*2);
     else
         subplot_ax_spikes = [];
         subplot_ax_speed=[];
     end
    plotExampleCellShift(data,cluID,ops,subplot_ax_spikes,subplot_ax_speed);

    
end
