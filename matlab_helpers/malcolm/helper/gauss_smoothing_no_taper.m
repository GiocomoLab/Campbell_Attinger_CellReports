function x=gauss_smoothing_no_taper(x,smoothWindow,smoothSigma)
% MGC 10/4/2017
if smoothWindow>0 && smoothSigma>0
    % define gaussian filter
    gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
    % convolve with data
    x = conv([repmat(x(1),floor(smoothWindow/2),1); x; repmat(x(end),floor(smoothWindow/2),1)],gauss_filter,'valid');
end