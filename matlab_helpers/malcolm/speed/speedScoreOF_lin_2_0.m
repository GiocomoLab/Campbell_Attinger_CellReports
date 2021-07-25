function [speedScore,params,hh] = speedScoreOF_lin_2_0(session,tetrode,cluster,boxSize)
% function to calculate speed score in OF
% based on Kropff et al 2015
% Malcolm Campbell, 4/1/16

% threshold for speed
speedThreshold = 2;
smoothSpeedTrace = 1;
makePlots = 1;

% load the position and spike files
cell = regexprep(strcat('T',num2str(tetrode),'C',num2str(cluster)),'[^\w'']','');
spikefile = strcat(session,'_',cell,'.mat');
posfile = strcat(session,'_pos.mat');
load(spikefile);
load(posfile);

% take out NaN's and replace them with neighboring values
positions = {posx, posy};
for k = 1:2
    pos_temp = positions{k};
    nan_ind = find(isnan(pos_temp));
    for m = 1:numel(nan_ind)
        if nan_ind(m) - 1 == 0
            temp = find(~isnan(pos_temp),1,'first');
            pos_temp(nan_ind(m)) = pos_temp(temp);
        else
            pos_temp(nan_ind(m)) = pos_temp(nan_ind(m)-1);
        end
    end
    positions{k} = pos_temp;
end
posx = positions{1}; posy = positions{2};

% Scale the coordinates using the shape information
minX = nanmin(posx); maxX = nanmax(posx);
minY = nanmin(posy); maxY = nanmax(posy);
xLength = maxX - minX; yLength = maxY - minY; sLength = max([xLength, yLength]);
scale = boxSize / sLength;
posx = posx * scale; posy = posy * scale;

% compute speed at every time point
velx = diff([posx(1); posx]); vely = diff([posy(1); posy]); dt = 0.02;
speed = sqrt(velx.^2+vely.^2)/dt;
speed(speed > 100) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
if smoothSpeedTrace
    speed = gauss_smoothing(speed,10); % gaussian kernel smoothing
end

% estimate instantaneous firing rate by smoothing spike count histogram
h = hist(cellTS,post);
fr = gauss_smoothing(h,20);  % smooth with gaussian kernel (sigma=20)
select = speed > speedThreshold; % throw out speeds below threshold

% calculate speed score by correlating instantaneous firing rate and speed
speed = reshape(speed,numel(post),1);
fr = reshape(fr,numel(post),1);
speedFilt = speed(select);
frFilt = fr(select)*50; % multiply by 50 to change to units of Hz
speedScore =corr(speedFilt,frFilt);
params = [ones(sum(select),1) speedFilt]\frFilt;

if makePlots
    edges = 0:5:max(speedFilt);
    [~,bin] = histc(speedFilt,0:5:max(speedFilt));
    frMeans = nan(max(bin),1);
    frErrs = nan(max(bin),1);
    for i = 1:max(bin)
        frMeans(i) = mean(frFilt(bin==i));
        frErrs(i) = std(frFilt(bin==i)); %/sqrt(sum(bin==i));
    end
    hh = figure('Visible','off'); shadedErrorBar(edges(1:end-1)+2.5,frMeans,frErrs,'k');
    title(sprintf('OF: speed score = %0.2f, slope = %0.2f',speedScore,params(2)));
    xlim([0 max(speedFilt)])
    ylim([0 max(frMeans)+max(frErrs)+1])
    h = refline(params(2),params(1));
    h.LineStyle='--';
    h.LineWidth=2;   
else
    hh = figure;
end

end