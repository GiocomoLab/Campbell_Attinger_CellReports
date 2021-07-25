%% run this first
loop_tbt_decode_position
%%
matfiles = dir('loop_tbt_output_dir');
errorMat = [];
timeError = [];
region = {};
x_vec = 1:2:399;
ops = load_default_opt;
ops.BinWidth = 2;
ops.edges = 0:ops.BinWidth:400;
ops.xbinedges = ops.edges;
ops.xbincent = .5*ops.edges(1:end-1)+.5*ops.edges(2:end);
savefig  = false;


for iF=1:numel(matfiles)
    [~,sn]=fileparts(matfiles(iF).name);
    
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    errorMat = cat(3,errorMat,squeeze(data_out.scoreMat2));
    
    t1=data_out.time_error*-1;
    
    gain_trial = ismember(data_out.trials,data_out.trials(1)+[7:10]);
    timeError = cat(1,timeError,nanmean(t1(gain_trial)));
    
    a=data_out.region;
    region=cat(1,region,a(1));
    
end

%%

%%
regs = {'MEC','VISp','RS'};
ALL={};
ALL_ANOVA=[];
for ii=1:3
    idx = startsWith(region,regs{ii});
    tmpBL = squeeze(nanmean(nanmean(abs(errorMat(1:6,:,idx)),1),2));
    tmpGC = squeeze(nanmean(nanmean(abs(errorMat(7:10,:,idx)),1),2));
    tmp = squeeze(tmpGC)-squeeze(tmpBL);
    ALL{ii}=tmp;
    ALL_ANOVA = cat(1,ALL_ANOVA,[tmp,ones(size(tmp))*ii]);
end

figure
plotSpread(ALL)
xticklabels(regs)
ylabel(' |Error| difference')
hold on
boxplot(ALL_ANOVA(:,1),ALL_ANOVA(:,2),'PlotStyle','compact')
xticklabels(regs)
%%
%%
regs = {'MEC','VISp','RS'};
ALL={};
figure
for ii=1:3
    subplot(1,3,ii)
    idx = startsWith(region,regs{ii});
    tmpBL = squeeze(nanmean(nanmean(abs(errorMat(1:6,:,idx)),1),2));
    tmpGC = squeeze(nanmean(nanmean(abs(errorMat(7:10,:,idx)),1),2));
    scatter(tmpBL,tmpGC,155,'k.')
    axis image
    xlim([0 110])
    ylim([0 110])
    xlabel('Baseline error')
    ylabel('Gain error')
    yticks(0:50:100)
    xticks(0:50:100)
    title(regs{ii})
    
end


%%
