function [y] = preprocess(x,win_width,perc_overlap, tresh,p)
%preprocess.m is a function that applies de-cklicking by applying a  
%weighting function and median filtering to the signal. 
% This processing is applied within each sliding window.
% The function returns the de-clicked time domain signal.
%
% INPUTS:
% - x = time domain signal (each channel is in a separate column)
% - win_width = window length in samples (processing is carried out within
%               each window)
% - pec_overlap = percentage overlap between consecutive window (e.g. 0.5
%               is 50% overlap)
% - tresh = threshold used in weighting function 
%           (see Gillespie et al (2013), Section II.A.1 for details)
% - p = power factor used in weighting function
%           (see Gillespie et al (2013), Section II.A.1 for details)
%
% OUTPUTS:
% y = de-clicked time domain signal

% Written by Pina Gruden, RCUH, Hawaii, 2021.


% References:
% -  D. Gillespie, M. Caillat, J. Gordon, and P. R. White (2013). Automatic 
% detection and classication of odontocete whistles, J. Acoust. Soc. Am. 
% 134(3), 2427-2437.
% - P. Gruden, E-M Nosal, and E. Oleson (2021 In revision). Tracking Time 
% Differences of Arrivals of multiple sound sources in the presence of 
% clutter and missed detections, J. Acoust. Soc. Am.


slide_incr=round(win_width*(1-perc_overlap)); % window advancment (increment)
w=triang(win_width);w=w(:);

numstps = ceil((length(x)-win_width)/slide_incr); %Number of windows 
start=1;
y = zeros(size(x));

for n=1:numstps %do sliding window
    xi= x(start:start+win_width-1,:);
    
    %----------------Remove clicks from the signal---------------------
    %Section II.A.1 in Gillespie et al. (2013)
    sd=std(xi);
    mu=mean(xi);
    wi=1./(1+((xi-mu)/(tresh*sd)).^p);% Weighting function
    xn = xi.*wi;  % weighted signal (declicked)

    %-----------FFT, median filter and IFFT-----------
    Y = fft(xn.*w,win_width); % FFT with hanning window
    L=length(Y);
    
    Ymag= 20*log10(abs(Y(1:(L/2)+1,:))); % Magnitude squared on a dB scale, single sided including Nyquist bin
    Ymedianfilt_dB = medfilt1(Ymag,61);
    Ymedianfilt = 10.^(Ymedianfilt_dB/20); %magnitude on a linear scale
    Yhat_single = Y(1:(L/2)+1,:)./Ymedianfilt; %single sided
    Yhat = [Yhat_single; conj(flipud(Yhat_single(2:(L/2),:)))];
 
    y(start:start+win_width-1,:)=y(start:start+win_width-1,:)+real(ifft(Yhat,[],1));
    
    
    start=start+slide_incr;
end


end

