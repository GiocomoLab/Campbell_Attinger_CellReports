close all;
n_mice = 1;

show_lm_bar = 1;
show_attract_ball = 1;

lmark_strength_list =2./15;

delta_x = (2* pi)/120. * (6);


noise_list = [0 0.0095*ones(1,26)]; % found that this noise value produces ok plots

attract_k_list = fliplr((2 * pi/40) * [1.0])*ones(1,1);


lmark_k_list = (2 * pi/40) + ( attract_k_list.* [1.15]);

bl_k = 2*pi/40;
lmark_hetero_list = 0;

phaseplot=figure();
shiftplot=figure();
x_list  = 0:delta_x:600;

landmark_phase = 2*pi;
attract_phases = 2*pi;
attract_phases_bl = 2*pi;
attract_bl = zeros(1,length(x_list));
save_video = false;
plot_shift = false;

% for xcorr
n_reps = length(noise_list);
stride = 40;
width = 100;
onsets = 1:50:(length(x_list)-width);
centers = x_list(onsets);
centers = centers+x_list(width)/2;
n_steps = numel(onsets);
PEAKS = zeros(n_reps,n_steps);
LAGS=PEAKS;
cmap = cbrewer('qual','Set1',numel(noise_list));
shift_fig = figure();

%compute 'baseline map', i.e. attractor phase if there would have been no
%gain change
for iT = 1:length(x_list)
    x = x_list(iT);
    %comment second part here for noise free version
    landmark_phase = landmark_phase+ delta_x * lmark_k_list(1);% + sqrt(delta_x * 0.001)*randn(1);
    attract_phases_bl = attract_phases_bl + delta_x * attract_k_list(1);
    attract_phases_bl = attract_phases_bl + delta_x * sin(landmark_phase - attract_phases_bl)*lmark_strength_list(1);
    attract_bl(iT)=landmark_phase;
    
    
    mean_attract_baseline(iT) = (mean(exp(j*attract_phases_bl)));
    mean_landmark_baseline(iT) = exp(j*landmark_phase);
end


lmark_hetero=0;
noise=0;
gain_onset = 100;
for i_lm_strength = 1:length(lmark_strength_list)
    lmark_base_strength = lmark_strength_list(i_lm_strength);
    for i_K = 1:length(attract_k_list)
        MEAN_SHIFT = zeros(length(noise_list),length(x_list));
        for i_noise = 1:length(noise_list)
            lmark_base_strength = lmark_strength_list(i_lm_strength)+(rand(1)-0.5)*0.1;
            
            attract_k = attract_k_list(i_K);%+randn(1)*0.5;
            landmark_k = lmark_k_list(i_K);
            
            attract_phases = 0 * sign(randn(n_mice, 1))   + 2 * pi;
            
            landmark_phase = 0 + 2*pi;
            attract_noise = noise_list(i_noise);
            
            
            reset_val = 1;
            
            for iT = 1:length(x_list)
                x = x_list(iT);
                
                landmark_phase = landmark_phase + delta_x * landmark_k;
                if x_list(iT)<gain_onset
                    attract_phases = attract_phases+delta_x*landmark_k;
                else
                    attract_phases = attract_phases + delta_x * attract_k;
                end
                cur_landmark_strength = lmark_base_strength * (1+ lmark_hetero * cos(landmark_phase));
                
                attract_phases = attract_phases + delta_x * sin(landmark_phase - attract_phases)*lmark_base_strength*reset_val;
                %noise_free = attract_phases;
                attract_phases = attract_phases + sqrt(delta_x * attract_noise) * randn(n_mice, 1);
                
                
                dd = imag(log(exp(j*attract_phases) .* conj(exp(j*landmark_phase))));
                % threshold for 'remapping', if shift hits this,
                thresh  = -pi*.5; %
                if dd<thresh% 
                    if i_noise > 1
                        noise = randn(1);
                    else
                        noise =0;
                    end
                    attract_phases =  2 * pi+noise;
                    %
                    landmark_phase = 2*pi+noise;
                    landmark_k = bl_k;
                    %
                    attract_k = bl_k;
                    reset_val = 0;
                end
             
                
                coherence(iT) = abs(mean(exp(1i*attract_phases)))^2; %Average dot product between two different firing rates
                mean_attract_phase(iT) = (mean(exp(1i*attract_phases)));
                mean_landmark_phase(iT) = exp(1i*landmark_phase);
                
            end
            mean_coherence = mean(coherence(x_list>max(x_list/2)));
            mean_shift_vs_t = imag(log(mean_attract_phase .* conj(mean_landmark_phase)));
            mean_shift = mean(mean_shift_vs_t(x_list>max(x_list/2)));
            mean_shift_vs_baseline = imag(log(mean_attract_phase .* conj(mean_landmark_baseline)));
            MEAN_SHIFT(i_noise,:)=mean_shift_vs_baseline;
            %MEAN_SHIFT(i_noise,:)=mean_shift_vs_t;
            fprintf(' AttractK %f, noise %f, Hetero %f  has ||| shift %f, coherence %f,  \n', attract_k, attract_noise, lmark_hetero, mean_shift, mean_coherence);
            %close all;
            mean_coherence_square(i_K, i_noise) = mean_coherence;
            mean_shift_square(i_K, i_noise) = mean_shift;
            
            
  
            
            
            
            figure(phaseplot)
            subplot(2,1,1)
            hold on
            plot(x_list, .5 * (1+ real(mean_attract_phase)), 'r', 'linewidth', 5);
            plot(x_list,.5* (1+real(mean_landmark_phase)),'b--','linewidth',3)
            title('mean phase')
            legend({'Attractor','Landmark'})
            %ylim([-1 3]);
            xlim([-5, 35]);
            subplot(2,1,2)
            hold on
         
            [ff,lags]=xcorr(mean_attract_phase(1100:end),mean_landmark_baseline(1100:end),'coeff');
            hold on
            plot(lags,real(ff))
            xlim([-20 20])
            for iStep =1:n_steps
                idx = onsets(iStep)-1+(1:width);
                [ff,lags]=xcorr(mean_attract_phase(idx),mean_landmark_baseline(idx),'coeff',20);
                [ma,mi]=max(real(ff));
                PEAKS(i_noise,iStep)=ma;
                LAGS(i_noise,iStep)=lags(mi);
            end
            
        end
    end
end

%%

%onsets = arrayfun(@(ROWIDX) strfind(MEAN_SHIFT(ROWIDX,:)<thresh,[0 1]), 1:size(MEAN_SHIFT,1),'UniformOutput',false);
remaps = [];
remap_times = [];
noise_free = [];
take_idx = -200:200;
for iR=1:size(MEAN_SHIFT,1)
    onsets = strfind(diff(MEAN_SHIFT(iR,:))>pi/8,[0 1]);
    for iO=1:min(1,numel(onsets))
        if onsets(1)+take_idx(end)>size(MEAN_SHIFT,2)
            continue
        end
        remap_times(end+1) = gain_onset-x_list(onsets(1));
        remaps = cat(1,remaps,MEAN_SHIFT(iR,onsets(1)+take_idx));
        noise_free = cat(1,noise_free,MEAN_SHIFT(1,onsets(1)+take_idx));
    end
end
figure('Renderer','Painters','Position',[680   421   248   435])
subplot(2,1,1)
mm=-1*median(remaps);
ss=-1*mad(remaps,1);
boundedline(take_idx*delta_x,mm,ss,'alpha')
hold on
mm=-1*mean(noise_free);
ss=std(noise_free)/sqrt(size(remaps,1));
%boundedline(take_idx/delta_x,mm,ss,'alpha','cmap',[.6 .6 .6])
%yline(noise_free_level)
axis tight
ylabel('Attractor shift')
xlabel('delta x')
xl=[-50 0.1]
ylim([-pi/8 pi])
xlim(xl)
box off
subplot(2,1,2)
plotSpread(remap_times','xyOri','flipped','spreadWidth',.1,'spreadFcn',{'xp',.1})
boxplot(remap_times,'PlotStyle','compact','Orientation','horizontal')

xlim(xl)
box off