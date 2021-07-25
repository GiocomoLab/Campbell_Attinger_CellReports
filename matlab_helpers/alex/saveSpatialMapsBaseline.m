%% setting parameters
ops = load_default_opt;
ops.xbin = 2;
ops.xbinedges = 0:ops.xbin:400;
ops.xbincent = ops.xbinedges(1:end-1)+ops.xbin/2;
ops.nBins = numel(ops.xbincent);
lw=ceil(3*ops.smoothSigma_dist);
wx=-lw:lw;
gw=exp(-wx.^2/(2*ops.smoothSigma_dist^2)); gw=gw/sum(gw);
ops.edges = ops.xbinedges;
ops.filter = gw';
ops.trials=5:20;
ops.num_pcs_decoding = 10;
regions = {'MEC','VISp','RS'};

%% setting paths and getting the session files
paths = get_paths;
data_path = paths.data;
savepath = fullfile(paths.data,'spatialMaps');

if ~isfolder(savepath)
    mkdir(savepath)
end
files = dir(fullfile(data_path,'*.mat'));
animal_list = {};
valid_files = {};
for iF=1:numel(files)
    [~,sn]=fileparts(files(iF).name);
    parts = strsplit(sn,'_');
    animal = parts{1};
    date = parts{2};
    animal_date = [animal, '_',date];
    if ~contains(sn,{'mismatch','playback','dark'}) && ~ismember(animal_date,animal_list) %use all files except mismatch, playback, dark
        valid_files{end+1}=files(iF).name;
        animal_list{end+1}=animal_date;
    end
end
%%
% p=gcp('nocreate');
% if isempty(p)
%     parpool(12);
% end

%% loop to save spatial maps
for iF=1:numel(valid_files)
    data = load(fullfile(data_path,valid_files{iF}));
    [~,sn]=fileparts(valid_files{iF});
    if ~isfield(data,'anatomy')
        continue
    end
    if isfield(data.anatomy,'parent_shifted')
        reg = data.anatomy.parent_shifted;
    else
        reg = data.anatomy.cluster_parent;
    end
    if iscolumn(reg)
        reg=reg';
    end
    %only use good_units that are part of MEC, VIS or
    good_cells_idx = data.sp.cgs == 2 & startsWith(reg,regions);
    if nnz(good_cells_idx)<2
        continue
    end
    if ~all(ismember(ops.trials,data.trial))
        continue
    end
    
    tmp_reg = reg(good_cells_idx);
    good_cells = data.sp.cids(good_cells_idx);
    good_idx = ismember(data.sp.clu,good_cells);
    clu_tmp = data.sp.clu(good_idx);
    st_tmp = data.sp.st(good_idx);
    
    [a,~,clus]=unique(clu_tmp);
    if isfield(data.anatomy,'depth')
        depth = data.anatomy.depth(good_cells_idx);
    else
        if ~isnan(data.anatomy.z2)
            depth = data.anatomy.tip_distance(good_cells_idx)-data.anatomy.z2;
        else
            depth = data.anatomy.tip_distance(good_cells_idx)-1000;
        end
    end
    if ~isrow(depth)
        depth=depth';
    end
    
    nClu = numel(a);
    ops_temp = ops;
    ops_temp.trials = 1:max(data.trial);
    [corrMat,frMat,~]=trialCorrMat(good_cells,ops.trials,data,ops);

   

   
    
    
    
    
    mf = matfile(fullfile(savepath,sn),'Writable',true);
    mf.cluID = good_cells;
    mf.region = tmp_reg;
    mf.depth = depth;
    mf.frMat = frMat;
    mf.corrMat = corrMat;
end