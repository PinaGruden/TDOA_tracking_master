function [RT,lags,RT_envelope,t] = gcc(x,y,fs,win_width,perc_overlap,Wf, freq_filter, HT)
%gcc.m function computes the generalized cross-correlation (GCC) of two signals. 
%Three types of GCC are available: standrad, phase transform, and smoothed
%coherence (see paper by Carter 1987). In addition the function returns an
%envelope of the cross-corerlation computed through the Hilbert transform.
%
% INPUTS:
% - x = data- each column is a channel, at the moment just two channels are
% considered, so x should be a M x N matrix where M is the samples and N is
% the channel
% - fs = sampling frequency
% - win_width = length of the window in samples. This MUST be even number.
% - perc_overlap = percentage of overlap between windows (e.g. for 75%
% overlap between windows input 0.75)
% - Wf - type of frequency weighting method (see Carter (1987))- choose between:
      %~ 'scc' - standard cross correlation
      %~ 'phat' - phase only transform (it whitenes the cross-spectrum)
      % 'scot' - smoothed coherence transform
% - freq_filter - Specified frequencies for band-pass filter. Default is
% not to perform filtering.
% - HT - logical - 1 for computing Hilbert transform and returning the
% envelope of the cross-correlation, 0 for not returning the envelope.
% Default is not to compute the envelope.
%
% OUTPUTS:
% - RT = cross-correlogram - each column is a
% cross-correlation of the two signals for a given window. It is a M x N
% matrix, where M depends on the window length (win_width) and N depends on
% the number of frames (numsteps)
% - lags = time lags 
% - RT_envelope = A matrix containing envelopes of the cross-correlations.
% Each column is an envelope for the cross-correlation in a given window.
% It is a M x N matrix,  where M depends on the window length (win_width) and N depends on
% the number of frames (numsteps)
% - t = vector of start times for each frame
%
% Given the outputs RT and lags, and assuming one source
% the time delay tau between channels would be computed as: 
%[~,idx] = max(abs(RT),[],1);
% tau = lags(idx);
%
%I am using definition for cross correlation:
%R_{xy}(\tau) = E[{x}(t){y}(t-\tau)]
%Thus the Cross-spectrum is:
%S_{xy}(f) = \lim_{T \rightarrow \infty} \frac{E[X_T(f)Y_T^*(f)] }{T}
% this means that if I get a NEGATIVE \tau the Y signal IS DELAYED
% with respect to X.
%!!!!!!!!!

%Pina Gruden, July 2019


if nargin<8, HT = 0; end
if nargin<7,freq_filter=0; HT = 0; end

assert(~mod(win_width,2),'The window width must be even')

%~~~~~~~~~~ Detrend the data ~~~~~~~~~~~~~~~~~~~~~~
x=detrend(x);
if ~isempty(y)
    y=detrend(y);
end
%~~~~~~~~~SPECIFY SLIDING WINDOW PARAMS~~~~~~~~~~~~~

slide_incr=round(win_width*(1-perc_overlap)); % window advancment (increment)- in samples
w=hann(win_width);w=w(:);
% w=chebwin(win_width,60);w=w(:); %w = chebwin(L,r) returns an L-point Chebyshev window using sidelobe magnitude factor r dB.

dt = slide_incr/fs; %time increment between frames in seconds
t=0:dt:(length(x)/fs-dt); %start times of frames, in seconds (for easier plotting)

% numstps = ceil((length(x)-(win_width))/slide_incr); %Number of windows
%If I wanted to include the last few samples/the whole of x, I would need to do:
%numstps = ceil((length(x)-(win_width))/slide_incr) +1; which increases the
%numsteps by 1 (this is due to Matlab starting indexing with 1 and not 0),
%then I would also have to zeropad the x -- so:
% numstps = ceil((length(x)-(win_width))/slide_incr) +1; %WRONG
numstps = ceil(length(x)/slide_incr);
nzeros=((numstps-1)*slide_incr+win_width)-length(x); %number of samples that is missing at the end to have a full frame (win_width long)
if isempty(y) %if there is no other file after x (i.e. end of encounter)
    x= [x;zeros(nzeros,size(x,2))]; %zeropad
else
    x= [x;y(1:nzeros,:)];%Add the required number of samples from y to the current file x
end


%~~~~~~~~~~~~~~~~DO A SLIDING WINDOW and COMPUTE GCC-PHAT~~~~~~~~~~~~~~~~

%preallocate and reset variables
start=1;
cols=[1,2]; %which channels to analyse
% nfft=win_width;
N=2*win_width-1; %lenght of cross-correlation
nfft=2^nextpow2(N); %before taking the FFT one needs to zero pad the signal
%in order to correct for the warap-around (circular convolution) effects- the cross-correlation is
%expected to return 2*N-1 length vector, but transposing to freq domain
%results in N length vector, so zero pad to correct length.
RT=zeros(nfft,numstps);
RT_envelope=RT;

% figure
for k=1:numstps %do sliding window
    xi= x(start:start+win_width-1,cols);
    
    %     %////////////// De-click/////////////////////
    %     p=6; tresh=5;
    %     sd=std(xi);
    %     mu=mean(xi);
    %     wi=1./(1+((xi-mu)/(tresh*sd)).^p);% Weighting function
    %     xn = xi.*wi;  % weighted signal (declicked)
    %
    % %     subplot(211),plot(xi(:,1)), title('Original signal')
    % %     subplot(212),plot(xn(:,1)), title('De-clicked signal')
    % %     drawnow
    % %     pause
    %
    %     xi=xn;
    %     %//////////////////////////////////////////
    
    %-----------FFT-----------
    Y = fft(xi.*w,nfft,1); % FFT for each column for a given window
    
    
    % % %     %////////Median filter aross frequency//////////
    % %it is redundant with PHAT and SCOT since it's a whitening of sorts
    %     Ymag= 20*log10(abs(Y)); % Magnitude squared on a dB scale
    %     Ymedianfilt_dB = medfilt1(Ymag,61);
    %     Ymedianfilt = 10.^(Ymedianfilt_dB/20); %magnitude on a linear scale
    %     Yhat = Y./Ymedianfilt;
    %     Y=Yhat;
    % %This does not affect SCOT or PHAT, but it makes SCC perform similar to
    % %PHAT and SCOT
    %     %//////////////////
    if freq_filter~=0
        %--------------- Set frequency range -----------------
        f0=freq_filter(1);
        f1=freq_filter(2);
        df=fs/nfft;
        freq = 0:df:fs-df;
        indx=find((freq > f0) & (freq<f1) | (freq > freq(end)-f1) & (freq<freq(end)-f0));
        
        %         range= (freq > f0) & (freq<f1) | (freq > freq(end)-f1) & (freq<freq(end)-f0);
        %         Sxy=Sxy.*range';
        % %         Sxx=Sxx.*range';
        % %         Syy=Syy.*range';
    else % no frequency filtering
        indx=[1:length(Y)];
    end
    
    %---------Cross-Spectrum Sxy-------------
    Sxy_temp=Y(indx,1).*conj(Y(indx,2)); % I could correct it for the duration or obtain
    %unbiased estimate, but I dont think it's necessary for this application
    
    %NOTE: since I'm complex conjugating the second channel-
    % when tau is computed if it is negative it means that signal in 2nd
    % channel is delayed with respect to 1 channel.
    
    %compute auto-spectra
    Sxx=Y(indx,1).*conj(Y(indx,1));
    Syy=Y(indx,2).*conj(Y(indx,2));
    %---------------GCC-----------------------
    switch Wf
        case 'scc'
            %            Sxy=Sxy; %stays the same
        case 'phat'
            Sxy_temp=exp(1i*angle(Sxy_temp));
        case 'scot'
            Sxy_temp=Sxy_temp./sqrt(Sxx.*Syy);
        otherwise
            disp('Not a valid method')
    end
    
    %     % NOTE: before I was using the code below, but in some cases where
    %     % Sxy was super small or 0 it returned NaN when computing W - thus changed to the above code
    %     switch Wf
    %         case 'scc'
    %             W=1;
    %         case 'phat'
    %             W=1./abs(Sxy);  % added very small number (+1e-16) to avoid division by 0
    %         case 'scot'
    %             Sxx=Y(:,1).*conj(Y(:,1));
    %             Syy=Y(:,2).*conj(Y(:,2));
    %             W=1./sqrt(Sxx.*Syy); % added very small number (+1e-16) to avoid division by 0
    %         otherwise
    %             disp('Not a valid method')
    %     end
    %     Sxy=W.*Sxy;
    %
    %     %NOTE: PHAT can also be computed as Sxy=exp(1i*angle(Sxy)); which is
    %     %the way it's computed in Matlab's inbuilt gccphat. It is the same as
    %     % the above ((1/abs(Sxy))*Sxy) since Sxy=abs(Sxy)*exp(j*phi) and
    %     %Swhite = W*Sxy = (1/abs(Sxy))*abs(Sxy)*exp(j*phi)=exp(j*phi)
    
    
    % Now set the vector Sxy to be the same length as Y, and with the values
    % outside desired frequency range to be zero:
    Sxy =zeros(size(Y,1),1);
    Sxy(indx)=Sxy_temp;
    
    
    %--------IFFT to obtain cross-correlation Rxy(tau)--------
    
    sxy= 2*ifft(Sxy(1:nfft/2),nfft); % this gives a complex signal (need for Hilbert transform)
    Rxy=fftshift(real(sxy),1);  %Eq. (13.84) from Bendant and Piersol (2010), Random data
    
    %     Rxx=ifft(Sxx,nfft);
    %     Ryy=ifft(Syy,nfft);
    %     Rxx0=Rxx(1);
    %     Ryy0=Ryy(1);
    %     sigx=sqrt(Rxx0);
    %     sigy=sqrt(Ryy0);
    
    % use the ifftshif: the Rxy consist of
    %values for positive lags, followed by values for negative
    %lags- by using fftshift, you swap halves of Rxy (each column),
    %so that the values are first for negative, then positive lags.
    % Due to numerical errors, the ifft(Sxy) returns a complex
    % number instead of real, so take only the real part real(ifft(Sxy))
    
    %-------------Compute the envelope of the Rxy(tau)-----------
    if HT == 1
        Rxy_tilde =fftshift(imag(sxy),1); %Eq. (13.84) from Bendant and Piersol (2010), Random data
        Rxy_envelope = sqrt(Rxy.^2 + Rxy_tilde.^2); % Eq. (13.73) from Bendant and Piersol (2010), Random data
        RT_envelope(:,k) = Rxy_envelope;
        %         RT_envelope(:,k) = Rxy_envelope./(sigx*sigy); %normalized
    end
    
    %-------Save Rxy(tau) for a given time step k-------
    %this builds a "cross-correlogram"
    RT(:,k)= Rxy;
    %     RT(:,k)= Rxy./(sigx*sigy); %normalized- correlation coefficient
    
%             figure
%             subplot(411)
%             plot(xi(:,1),'k'), title('1st channel')
%             xlabel('Samples'),ylabel('Amplitude')
%             subplot(412)
%             plot(xi(:,2),'k'), title('2nd channel')
%             xlabel('Samples'),ylabel('Amplitude')
%             subplot(413)
%             lags1=(-nfft/2:nfft/2-1)/fs;
%             plot(lags1,Rxy),title('Cross-correlation')
%             xlabel('Lags (s)'), ylabel('R_{xy}')
%             subplot(414)
%             plot(lags1,Rxy_envelope),title('Envelope of Cross-correlation')
%             xlabel('Lags (s)'), ylabel('A_{xy}')
%     %         xlim([-0.002,0.002])
%             drawnow
%             pause
    %
    start=start+slide_incr;
end

% truncate the RT to be win_width long, centered at 0 lag:
RT = RT((nfft/2+1)-win_width/2:(nfft/2+1)+win_width/2-1,:);
if HT==1
    RT_envelope = RT_envelope((nfft/2+1)-win_width/2:(nfft/2+1)+win_width/2-1,:);
end
%Define lags
lags= -win_width/2:win_width/2-1; %Taking only win_width points in the middle of cross-correlation
lags=lags/fs;

end

