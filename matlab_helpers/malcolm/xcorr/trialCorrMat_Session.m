function corrMat = trialCorrMat_Session(session_name,opt)

if ~exist('opt','var')
    opt = load_default_opt;
end

paths = get_paths(opt);
dat = load(fullfile(paths.data,session_name));
cell_info = load(fullfile(paths.cell_info,opt.cell_file));
cell_info = cell_info.cell_info;

good_cells = cell_info.CellID(strcmp(cell_info.Session,session_name) & ...
    strcmp(cell_info.BrainRegion,opt.brain_region));

corrMat = trialCorrMat(good_cells,1:max(dat.trial),dat,opt);

end