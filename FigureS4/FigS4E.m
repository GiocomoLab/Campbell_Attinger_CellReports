% Makes the following figures:
% Fig S4E
% To use: modify the paths in the get_paths function to point to the right
% place on your machine.
% MGC 7/19/2021

%restoredefaultpath;

% add helper functions to path
%addpath(genpath('..\matlab_helpers\malcolm'));

% get paths
paths = get_paths; % *** NOTE: modify the paths in this function as appropriate ***

% load cell_info table
load(fullfile(paths.intermediate_data,fullfile('cell_info','cell_info_Apr2020'))); 

% analysis options
opt = load_default_opt;

% cells to plot
cells = {{'npH4_0319_contrast_2',262},...
    {'npJ1_0521_contrast_1',710},...
    {'npJ2_0512_contrast_1',337},...
    {'npJ5_0505_contrast_1',396},...
    {'npI4_0421_contrast_2',581},...
    {'npI3_0421_contrast_1',481},...
    {'npI1_0414_contrast_2',112},...
    {'npH1_0313_contrast_2',140},...
    {'npH5_0324_contrast_1',231}};

trials = 1:230; % trials to plot

%% plot example cells in Figure S4E
for i = 1:numel(cells)
    cellidx = find(strcmp(cell_info.Session,cells{i}{1}) & cell_info.CellID==cells{i}{2});
    celltitle = sprintf('%s_c%d_%s_d=%d',cells{i}{1},cells{i}{2},...
        cell_info.BrainRegion{cellidx},round(cell_info.DepthFromSurface(cellidx)));
    dat = load(fullfile(paths.data,cells{i}{1}));
    spiket = dat.sp.st(dat.sp.clu==cells{i}{2});
    [~,~,spike_idx] = histcounts(spiket,dat.post);
    hfig = figure('Position',[200 200 300 700]); hold on;
    hfig.Name = celltitle;
    frMat = calcTrialFRMat(cells{i}{2},trials,dat,opt);
    colormap(cbrewer('seq','BuPu',20));
    cb = colorbar('Location','SouthOutside');
    ylabel(cb,'firing rate (Hz)'); 
    imagesc(frMat);
    title(celltitle,'Interpreter','none');
    ylim([min(trials) max(trials)]);
    yticks([]);
    xlim([0 200]);
    xticks([0 100 200]);
    xticklabels([0 200 400]);
    xlabel('cm');
    
    % patches indicating contrast
    for tr = trials
        patch([-15 0 0 -15],[tr tr tr+1 tr+1]-0.5,get_color(1,dat.trial_contrast(tr)),...
            'EdgeColor',get_color(1,dat.trial_contrast(tr)));
    end
    for tr = 20.5:10:220.5
        plot([0 200],[tr tr],'w-');
    end
    xlim([-15 200]);
end