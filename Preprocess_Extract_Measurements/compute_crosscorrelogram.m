function [Rxy_envelope_ALL, lags, t_serialdate,t]=compute_crosscorrelogram(folder1,folder2save2, parameters,param_signal)
% compute_crosscorrelogram.m is a function that computes a
% cross-correlogram between two channels in audio data.

%INPUTS:
% - folder1 - a string specifying path to the folder where audio data is
%           located
% - folder2save2 - a structure specifying paths to where data is stored to.
%                   Has two fields:
%                   ~'rawcrosscorr' (path to where cross-correlogram 
%                   will be stored);
%                   ~'finalresults' (path to where final results will
%                    be stored)
% - parameters -  a structure with at least 13 fields containing parameters  
%               for the encounter and processing. The following 10 fields
%               are required:
%       ~ parameters.encounter- a string specifying the encounter name (should
%         match the name in Array_info.csv)
%       ~ parameters.c - a scalar specifying the speed of sound
%       ~ parameters.channels - 1 x 2 vector specifying sensor channels to
%           be used for processing
%       ~ parameters.d - a scalar specifying sensor separation
%       ~ parameters.method - a string specifying method for generalized
%           cross-correlation (GCC) computation. Should be 'scc'/'phat'/'scot' 
%           (see gcc.m for more info)
%       ~ parameters.window_length_s - window length (in s) over which the 
%           cross-correlation is computed
%       ~ parameters.overlap - proportion of overlap between consecutive
%           windows for cross-correlogram computation (between 0 and 1)
%       ~ parameters.plotcrosscorr- a scalar to specify to plot the 
%           cross-correlogram or not (0=no,1=yes)
%       ~ parameters.saveworksp - a scalar to specify to save the 
%           cross-correlogram or not (0=no,1=yes)
%       ~ parameters.dt - time increment between consecutive GCC windows (in s)
% - param_signal - a structure with two fields:
%       ~ signal_type -  a string specifying the signal type ('whistles'/'clicks') 
%       ~ freq_filter - 1 x 2 vector specifying min and max frequency limits
%        for the signal_type (in Hz).

%OUTPUTS:
% - Rxy_envelope_ALL- cross-correlogram- a M x N matrix containing 
%   cross-correlation information, where M is number of TDOAs and N is
%   number of time steps.
% - lags - 1 x M vector of TDOAs,
% - t_serialdate - 1 x N vector of times (in a serial date format),
% - t - 1 x N vector of times (starting at 0 at the beginning of encounter).

%Pina Gruden, 2022, UH Manoa

signal_type =param_signal.signal_type;
freq_filter = param_signal.freq_filter;

channels=parameters.channels;
encounter=parameters.encounter;
window_length_s=parameters.window_length_s;
overlap=parameters.overlap;
dt=parameters.dt;
method = parameters.method;
d=parameters.d;
c=parameters.c;
plotfig=parameters.plotcrosscorr;
saveworksp=parameters.saveworksp;

file1 = dir(fullfile(folder1,'*.wav' ) ); % list all wav files in subfolder

%-----------------------------------------------------------------
%eliminate all '.' and '..' files:
%-----------------------------------------------------------------
indx=[];
for k = 1:size(file1,1) %iterate through files in the folder
    if file1(k).name(1) == '.'
        indx=[indx,k]; % skip '.' and '..' files
    end
end
file1(indx)=[];

N=size(file1,1);


%-----------------------------------------------------------------
% Get all date-time information 
%-----------------------------------------------------------------
%(this is important if there are missing times between consecutive files):
FormatIn='yyyyMMdd_HHmmss_SSS'; % use if youre using datetime function
datenumberString = getStrDateTime(file1,FormatIn);
timestamps = datetime(datenumberString,'InputFormat',FormatIn);

duration_file=zeros(N,1);
for k=1:N
info = audioinfo([file1(k).folder,'/',file1(k).name]);
duration_file(k)=info.Duration;
end

actual_gaps=diff(timestamps);
expected_gaps=seconds(round(duration_file(1:end-1))); 
% duration_file(1:end-1) - the last file not considered since there is nothing after it

index=actual_gaps>expected_gaps;
flagf=find(index)+1;
%need to add +1 since the index should be to the file before which the gap
%occurs.

%Number of missing time steps for each file:
Nstepsmissing=seconds(actual_gaps(index)-expected_gaps(index))./dt;

% gaps = [NaN; diff(timestamps)];
% 
% index=gaps>seconds(60.1);
% % if sum(index)>0
% %     fprintf('Warning:\n There are gaps in the recordings: before file %i.\n', find(index))  
% % end
% flagf=find(index); %this gives index of the file before which there is a gap in the recordings
% 
% % %Number of missing time steps for each file:
% % % dt=(1-overlap)*window_length_s;%time increment between frames in seconds
% Nstepsmissing=(seconds(gaps(index))-60)/dt;


%-----------------------------------------------------------------
% COMPUTE GCC
%-----------------------------------------------------------------
%pre-allocate
% between(timestamps(1),timestamps(end))
% M=(seconds(nansum(gaps))+60)/dt; %this assumes the last file is 60 s long.
% 
% [~,fs]=audioread([file1(1).folder,'/',file1(1).name]);

M=round((seconds(timestamps(end)-timestamps(1))+duration_file(end))./dt);

fs=info.SampleRate; %this assumes all files have the same sampling rate
L= round((d/c*1.1)*2*fs); %cross-correlogram lags are limited by the sensor 
% separation- the lags will be between -d/c and d/c and we add 10% on top 
% in case distance is not accurately measured

Rxy_envelope_ALL= zeros(L,M);

start=1;
    for k = 1:N %iterate through files in the folder
        
       
        %% Read the WAV file
        file=file1(k).name;
        disp(['Computing ', signal_type,' cross-correlogram for file ',file])
        [x,fs]=audioread([file1(k).folder,'/',file]);
        %Select the two channels to process
        x = x(:,channels);
        
        %read the next file (one window length of the next file)
        if k<N
        [x1,~]=audioread([folder1,'/',file1(k+1).name],[1,window_length_s*fs]); 
        x1 = x1(:,channels);
        else %you are reading the last file of the encounter, so there is no file after
          x1=[];  
        end
          
        %% Specify FILTERS and do FILTERING
        %~~~~~ Bandpass Butterworth ~~~~~~~~~
        [b,a] = butter(4,freq_filter/(fs/2),'bandpass'); %Coefficients for Butterworth
        x_filt = filtfilt(b,a,x); % Butterworth in forward-reverse fashion- so that it is zero phase
        x1_filt = filtfilt(b,a,x1);
        
        switch signal_type
            case 'whistles'
                %~~~~ Remove clicks and filter ~~~~~~~
                y = preprocess(x_filt,1024,0.5,5,6); %tested and it does not affect phase
                y1 = preprocess(x1_filt,1024,0.5,5,6);
            case 'clicks'
                y=x_filt;
                y1=x1_filt;
            otherwise
                disp('Not a valid signal type. Choose "whistles" or "clicks"')
        end

        
        %% SPECIFY PARAMS for GCC and do GCC
        
        win_length=window_length_s*fs; %WINDOW LENGTH FOR GCC in samples
        win_step=1-overlap; %percentage
        HT=1; %compute envelope via Hilbert transform
        
        [~,lags,RT_envelope,~] = gcc(y,y1,fs,win_length,overlap,method,freq_filter,HT);

        ind = find(abs(lags)<d/c*1.1); %d/c*1.1 add 10% more values to the search in case the separation wasn't measured accuretly
        lags=lags(ind);
        RT_envelope = RT_envelope(ind,:);

        [m,M]=size(RT_envelope);
        stop = start+M-1;
        Rxy_envelope_ALL(1:m,start:stop) = RT_envelope;
        if any(k==flagf-1) %if this is a file after which the gap occurs
            start = stop+1+round(Nstepsmissing(k==flagf-1));
        else
            start = stop+1;
        end

    end

    %trim away unused pre-allocated space
    Rxy_envelope_ALL=Rxy_envelope_ALL(1:m,:);

    % dt=((1-overlap)*win_length)/fs;%time increment between frames in seconds
    t=(0:(size(Rxy_envelope_ALL,2)-1)).*dt; %start times of frames, in seconds (for easier plotting)

    %determine duration of the last file:
%     info = audioinfo([file1(end).folder,'/',file1(end).name]);
%     dur_last=info.Duration;
    dur_last=duration_file(end);

%     t_datetime=timestamps(1):seconds(dt):timestamps(end)+seconds(dur_last)-seconds(dt);
    t_datetime=timestamps(1):seconds(dt):timestamps(end)+seconds(dur_last);
    t_serialdate=datenum(t_datetime);

    %% SAVE WORKSPACE containing CROSS-CORELOGRAM
    
    if saveworksp==1
    switch signal_type
        case 'whistles'
            %parameters_whistles=param_signal;
            C=whos('Rxy_envelope_ALL');
            if C.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
                save([folder2save2,encounter,'_whistles_rawCrossCorrelogram_ALL.mat'],...
                    'Rxy_envelope_ALL','lags','t_serialdate', 't',...
                    'parameters','param_signal','-v7.3')
            else
                save([folder2save2,encounter,'_whistles_rawCrossCorrelogram_ALL.mat'],...
                    'Rxy_envelope_ALL','lags','t_serialdate', 't',...
                    'parameters','param_signal')
            end
        case 'clicks'
            %parameters_clicks=param_signal;
            C=whos('Rxy_envelope_ALL');
            if C.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
                save([folder2save2,encounter,'_clicks_rawCrossCorrelogram_ALL.mat'],...
                    'Rxy_envelope_ALL','lags','t_serialdate', 't',...
                    'parameters','param_signal','-v7.3')
            else
                save([folder2save2,encounter,'_clicks_rawCrossCorrelogram_ALL.mat'],...
                    'Rxy_envelope_ALL','lags','t_serialdate', 't',...
                    'parameters','param_signal')
            end
        otherwise
            disp('Not a valid signal type. Choose "whistles" or "clicks"')
    end
    end

    if plotfig==1
        %f= figure('visible', 'off');
        figure;
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-d/c,d/c])
        xlabel('Local Time (HH:MM)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ',signal_type])
        set(gca,'FontSize',12)
        set(findall(gcf,'type','text'),'FontSize',14)
        
% %         caxis([0,0.1])
%%% Uncomment below if you want to save the cross-correlogram figure:
% f.PaperUnits='centimeters';
%         f.PaperPosition=[0,0,100,30];%size specified as panorama, 5:1- [left,bottom,width, height]
%         switch signal_type
%             case 'whistles'
%                 print([folder2save2,encounter,'_rawCrossCorrelogram_ALL_LocalTimeXaxis.jpg'],'-djpeg','-r0')
%             case 'clicks'
%                 print([folder2save2,encounter,'_rawCrossCorrelogram_ALL_LocalTimeXaxis.jpg'],'-djpeg','-r0')
%             otherwise
%                 disp('Not a valid signal type. Choose "whistles" or "clicks"')
%         end
%         close(f)
    end

end

