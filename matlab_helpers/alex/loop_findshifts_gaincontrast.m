

ops = load_default_opt;
ops.cm_ops = load_default_opt; %this one will get passed to the calcCorrMat (with 2cm spatial bin)
ops.factors = -.3:0.01:.3;

ops.SpatialBin = 1;
ops.idx = [10:ops.SpatialBin:390]/ops.SpatialBin;

ops.n_trials = 10;
ops.plotfig = false;

%% savedir =
savedir = fullfile('delay_estimate_gaincontrast');
imdir = fullfile(savedir,'images');
if ~isfolder(savedir)
    mkdir(savedir);
    
end

mf = matfile(fullfile(savedir,'parameters'),'Writable',true);
mf.ops = ops;





regions = {'VISp','MEC','RS'};
filenames = {};
fi = dir(fullfile(paths.data,'*gaincontrast*.m'));
for iF = 1:numel(fi)
    
filenames=cat(2,fullfile(fi(iF).folder,fi(iF).name));
end



%%
p=gcp('nocreate');
if isempty(p)
    parpool();
end

%%
zero_idx = find(ops.factors==0);

parfor iF=1:numel(filenames)
    [~,session_name]=fileparts(filenames{iF});
    
    if isfile(fullfile(savedir,[session_name '.mat']))
        disp('\t \t exisits')
        continue
    end
    try
        data = load(filenames{iF});
        ops_temp = ops;
        
        good_chunks = {};
        chunk_contrast = [];
        current_chunk = 1;
        chunk_size = [];
        tmp = [];
        current_contrast = data.trial_contrast(1);
        current_gain = data.trial_gain(1);
        start_new = false;
        for iT = 1:numel(data.trial_gain)
            bl_gain = data.trial_gain(iT) ==1;
            old_contrast = current_contrast;
            current_contrast = data.trial_contrast(iT);
            
            contrast_switch = current_contrast ~= old_contrast;
            
            old_gain = current_gain;
            current_gain= data.trial_gain(iT);
            gain_switch = old_gain ~= current_gain;
            
            
            if ~contrast_switch && bl_gain && ~gain_switch && numel(tmp)<ops.n_trials % add to list
                tmp(end+1)=iT;
            elseif contrast_switch %end of chunk
                %sprintf('contrast_switch at %d',iT)
                good_chunks{current_chunk}=tmp;
                chunk_size(current_chunk)=numel(tmp);
                chunk_contrast(current_chunk)=old_contrast;
                start_new = true;
            elseif gain_switch && ~bl_gain %gain change onset, save and reset tmp but don't start new
                %sprintf('gain_switch at %d',iT)
                
                good_chunks{current_chunk}=tmp;
                chunk_size(current_chunk)=numel(tmp);
                chunk_contrast(current_chunk)=current_contrast;
                tmp=[];
            elseif gain_switch && bl_gain %start new after gain comes back to baseline, no saving
                %sprintf('gain_switch at %d',iT)
                start_new = true;
            elseif numel(tmp)>=ops.n_trials
                %sprintf('max trials at %d',iT)
                good_chunks{current_chunk}=tmp;
                chunk_size(current_chunk)=numel(tmp);
                chunk_contrast(current_chunk)=old_contrast;
                start_new = true;
            end
            if start_new
                tmp=[];
                tmp(end+1)=iT;
                current_chunk = current_chunk+1;
                start_new = false;
            end
            
        end
        
        
        nC=nnz(data.sp.cgs==2);
        all_factors = nan(numel(good_chunks),nC);
        all_stability = all_factors;
        all_firingRate=all_factors;
        %ops_here.trial = find(data.trial_gain ==1 & data.trial_contrast==100);
        %run shift finding for each block
        
        for iRep=1:numel(good_chunks)
            ops_temp.trials = good_chunks{iRep};
            
            [data_out] = findBestShifts(data,ops_temp);
            [~,mi]=max(data_out.all_stability,[],2);
            factors = ops_temp.factors(mi);
            all_factors(iRep,:)=factors;
            all_stability(iRep,:)=data_out.stability;
            
        end
        
        
        mf =matfile(fullfile(savedir,session_name),'Writable',true);
        mf.chunk_contrast = chunk_contrast;
        mf.region = data_out.region;
        mf.subregion = data_out.sub_reg;
        mf.depth = data_out.depth;
        mf.CID = data_out.CID;
        mf.good_chunks = good_chunks;
        mf.all_factors = all_factors;
        mf.all_stability = all_stability;
        %mf.FR = all_firingRate;
    catch ME
        fprintf('%s \nFailed for %s: %d \n',ME.message,filenames{iF},iF)
        rethrow(ME)
    end
end
