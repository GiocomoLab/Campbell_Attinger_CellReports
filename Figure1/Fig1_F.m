paths= get_paths;
data_path = paths.data;
files = dir(fullfile(data_path,'*.mat'));
animal_list = {};
valid_files = {};
for iF=1:numel(files)
    [~,sn]=fileparts(files(iF).name);
    
    parts = strsplit(sn,'_');
    animal = parts{1};
    date = parts{2};
    animal_date = [animal, '_',date];
    if ~ismember(animal_date,animal_list)
        animal_list{end+1}=animal_date;
        valid_files{end+1}=files(iF).name;
    end
end

%%
N_UNITS = struct();
SPAN = struct();
for iF=1:numel(valid_files)
    data = load(fullfile(data_path,valid_files{iF}));
    if ~isfield(data,'anatomy')
        sprintf('no anatomy for %s \n',valid_files{iF})
        continue
    end
    if isfield(data.anatomy,'parent_shifted')
        region = data.anatomy.parent_shifted;
    else
        region = data.anatomy.cluster_parent;
    end
    
    [a,b,c]=unique(region);
    good_idx = data.sp.cgs'==2;
    for iR = 1:numel(a)
        this_idx = c==iR & good_idx;
        n_units = nnz(this_idx);
        [l,s]=bounds(data.anatomy.tip_distance(this_idx));
        dd = abs(diff([l,s]));
        if nnz(this_idx)>1 && isvarname(a{iR})
            if isfield(N_UNITS,a{iR})
                N_UNITS.(a{iR}) = cat(1,N_UNITS.(a{iR}),n_units);
                SPAN.(a{iR}) = cat(1,SPAN.(a{iR}),dd);
            else
                N_UNITS.(a{iR})=n_units;
                SPAN.(a{iR}) = dd;
            end
        end
        
    end
end
%%

fn = fieldnames(N_UNITS);
N_UNITS.RSC = [];
SPAN.RSC = [];
for ifn = 1:numel(fn)
    if startsWith(fn{ifn},'RS')
        N_UNITS.RSC = cat(1,N_UNITS.RSC,N_UNITS.(fn{ifn}));
        SPAN.RSC = cat(1,SPAN.RSC,SPAN.(fn{ifn}));
    end
end
%%
fn = fieldnames(N_UNITS);
for ifn = 1:numel(fn)
    idx = N_UNITS.(fn{ifn})<10;
    N_UNITS.(fn{ifn})(idx)=[];
    SPAN.(fn{ifn})(idx)=[];
end


%%

toplot = {'MEC','VISp','RSC'};
for iV =1:numel(toplot)
    tmp=N_UNITS.(toplot{iV});
    sprintf('%s: N= %.2f =- %.2f \n',toplot{iV},mean(tmp),std(tmp)/sqrt(nnz(tmp)))
end

%%
toplot = {'MEC','VISp','RSC'};
TMP = struct();
DENSITY = struct();
X=[];
G=[];
for iV =1:numel(toplot)
    TMP.(toplot{iV})=N_UNITS.(toplot{iV});
    tmp = N_UNITS.(toplot{iV})./SPAN.(toplot{iV});
    DENSITY.(toplot{iV}) = tmp;
    X = cat(1,X,tmp);
    G = cat(1,G,iV*ones(size(N_UNITS.(toplot{iV}))));
end
anova1(X,G)

toplot = {'MEC','VISp','RSC'};
TMP = struct();
DENSITY = struct();
X=[];
G=[];
for iV =1:numel(toplot)
    TMP.(toplot{iV})=N_UNITS.(toplot{iV});
    tmp = N_UNITS.(toplot{iV});
    DENSITY.(toplot{iV}) = tmp;
    X = cat(1,X,tmp);
    G = cat(1,G,iV*ones(size(N_UNITS.(toplot{iV}))));
end
anova1(X,G)
%%
cmap=cbrewer('qual','Set2',3,'pchip');
markers = {'o','x','v'};
toplot = {'MEC','VISp','RSC'};
TMP = struct();
DENSITY = struct();
X=[];
G=[];
figure('Position',[680   842   686   256])
for iV =1:numel(toplot)
    b_units=N_UNITS.(toplot{iV});
    span =SPAN.(toplot{iV});
    subplot(1,2,1)
    hold on
    scatter(b_units,span,25,cmap(iV,:),'o','filled')
    subplot(1,2,2)
    hold on
    scatter(b_units,span,25,[0 0 0],markers{iV})
end
subplot(1,2,1)
box off
legend(toplot)
legend('Location','southeast')
subplot(1,2,2)
box off
legend(toplot)
legend('Location','southeast')

