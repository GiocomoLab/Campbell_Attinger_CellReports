function [distScore,distPeak] = compute_dist_tuning(X,opt)
%COMPUTE_DIST_TUNING Computes a distance score using FFT
%   distScore is 3 x 1 vector: power of top three spectral components
%   distPeak is 3 x 1 vector: spatial period of these components
%   
%   X is a firing rate matrix, dimensions cells x distance bins

% Spatial Fourier transform
Fs = 1/opt.SpatialBin; % Sampling frequency                          
L = floor(size(X,2)/2)*2; % Length of signal (making sure it's even)
f = Fs*(0:(L/2))/L; % Define the frequency domain f
X = X(:,1:L);
X = (X-repmat(nanmean(X,2),1,size(X,2)))./repmat(nanstd(X,[],2),1,size(X,2)); % take z score
Y = fft(X');

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);

% min freq
min_freq_idx = find(f>opt.min_spatial_freq,1);
P1 = P1(min_freq_idx:end,:);
f = f(min_freq_idx:end);

% compute distance score
distScore = nan(size(X,1),3);
distPeak = nan(size(X,1),3);
for i = 1:size(X,1)
    [P1sort,sort_idx] = sort(P1(:,i),'descend');
    fsort = f(sort_idx);
    distScore(i,:) = P1sort(1:3);
    distPeak(i,:) = 1./fsort(1:3);
end

end