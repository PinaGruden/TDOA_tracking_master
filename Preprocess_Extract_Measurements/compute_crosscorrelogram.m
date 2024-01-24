function [Rxy_envelope_ALL, lags, t_serialdate,t]=compute_crosscorrelogram(folder1,folder2save2, parameters,param_signal)
% compute_crosscorrelogram.m is a function that computes a
% cross-correlogram between two channels in audio data.
%
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
%
%OUTPUTS:
% - Rxy_envelope_ALL- cross-correlogram- a M x N matrix containing 
%   cross-correlation information, where M is number of TDOAs and N is
%   number of time steps.
% - lags - 1 x M vector of TDOAs,
% - t_serialdate - 1 x N vector of times (in a serial date format),
% - t - 1 x N vector of times (starting at 0 at the beginning of encounter).
%
%
%Pina Gruden, 2022, UH Manoa

signal_type =param_signal.signal_type;
freq_filter = param_signal.freq_filter;

channels=parameters.channels;
encounter=parameters.encounter;
window_length_s=parameters.window_length_s;
overlap=parameters.overlap;
dt=parameters.dt;
method = parameters.method;
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
time_start = datetime(datenumberString,'InputFormat',FormatIn);

duration_file=zeros(N,1);
for k=1:N
info = audioinfo([file1(k).folder,'/',file1(k).name]);
duration_file(k)=info.Duration;
end
fs=info.SampleRate; %this assumes all files have the same sampling rate

% DETERMINE if there is a MISMATCH to expected length of files:
time_end= time_start + seconds(duration_file); 
%based on a duration of each file and the beginning time stamp this should
%be end of file timestamp and thus beginning of next file timestamp

starttiming_mismatch = zeros(N,1);
starttiming_mismatch(2:end)=seconds(time_start(2:end)-time_end(1:end-1));
% If it is a negative sign it means the file started that many milliseconds
% too early. If it is a positive sign it means the file started that many
% milliseconds/seconds too late and there is a gap


%-----------------------------------------------------------------
% Get GLOBAL TIME vector
%-----------------------------------------------------------------
find_lowerdt = @(x,y) floor(second(x)/y)*y;
%function find_lowerdt takes start file datetime information (as x) and
%time step information in seconds (as y).

find_upperdt = @(x,y) ceil(second(x)/y)*y;
%function find_upperdt takes start file datetime information (as x) and
%time step information in seconds (as y).

global_time_start_s= find_lowerdt(time_start(1),dt);%floor(second(time_start(1))/dt)*dt; %these are start seconds

newdatetime= @(x) datetime(year(x), month(x), day(x), hour(x), minute(x),0);

global_time_start= newdatetime(time_start(1))+ seconds(global_time_start_s);

last_file_end=time_start(end)+seconds(duration_file(end));

global_time_end_s= find_lowerdt(last_file_end,dt);%floor(second(last_file_end)/dt)*dt; % end seconds

global_time_end= newdatetime(last_file_end)+ seconds(global_time_end_s);

%start time of each time step from all files in datetime format
global_time=global_time_start:seconds(dt):global_time_end;

t_serialdate=datenum(global_time);


%-----------------------------------------------------------------
% COMPUTE GCC
%-----------------------------------------------------------------

tau_max_samples= round(parameters.tau_max*fs); %maximum physically possible lag in samples
window_length_samples=window_length_s*fs;

%Display a warning if max lag that is mathematically possible (i.e.
% dependent on the chosen window length) is shorter than
% lag that is physically possible
if window_length_s < parameters.tau_max
    fprintf(['WARNING: your chosen window length for cross-correlation ' ...
        'results in a maximum TDOA (parameters.window_length_s = %0.2f s) \n ' ...
        'that is shorter than what is physically possible (parameters.tau_max = %0.2f s).\n ' ...
        'Consider using a longer window in xcorrparam.win_width.\n'],parameters.window_length_s,parameters.tau_max);
end


%pre-allocate
L=min(tau_max_samples*2+1,window_length_samples*2-1);
M=length(global_time);

Rxy_envelope_ALL= zeros(L,M);

indx_timestep_start_k=zeros(N,1);
indx_timestep_end_k=zeros(N,1);

%Determine time steps where files start and end
for k = 1:N %iterate through files in the folder
    indx_timestep_start=find((global_time-time_start(k))<=0); %these are indices of all timesteps before the file starts
    indx_timestep_start_k(k)=indx_timestep_start(end); %take the last index in the list to be the time step where the file starts

    indx_timestep_end=find((global_time-time_end(k))<0); %these are indices of all timesteps before the file starts
    indx_timestep_end_k(k)=indx_timestep_end(end); %take the last index in the list to be the time step where the file starts
end


    for k = 1:N %iterate through files in the folder
      
        %% Read the WAV file
        file=file1(k).name;
        disp(['Computing ', signal_type,' cross-correlogram for file ',file])
        [x,fs]=audioread([file1(k).folder,'/',file]);
        %Select the two channels to process
        x = x(:,channels);

        % TIME CORECT FILES
        %Note- this assumes that files only overlap over one time step

        % Check for the beginning of file- if it starts at the same time
        % step that another file ends OR if it is alone (either first or
        % after a large gap

        % 1) Get x (input for GCC)

        [L1]=ismember(indx_timestep_start_k(k),indx_timestep_end_k);

        if L1 % the current file starts in the same time step when the previous file ended
               % the overlap is always sorted out at the end of the file,
               % so here we just chop off the beginning section that lies
               % int he overlap.

            % indx_overlap_timestep= indx_timestep_end_k(loc);
            start_time_read=global_time(indx_timestep_start_k(k)+1); %we want to start reading the file at the next time step (thus +1)
            start_sample=round(seconds(start_time_read-(time_start(k)))*fs)+1;

            x = x(start_sample:end,:);

           %indicate which time step it starts on:
            start_timestep=indx_timestep_start_k(k)+1;

        else % this file starts alone in its time step (so it is either a first file or file starting after a large gap ( >=dt or more))

            file_start_delay = seconds(time_start(k) - global_time(indx_timestep_start_k(k))); % delay in start of the file in seconds
            missing_samples=file_start_delay*fs; % delay in start of the file in samples
            x = [zeros(missing_samples,size(x,2));x];% zero pad the file
            % to start at the begininng of the time step

            %indicate which time step it starts on:
            start_timestep=indx_timestep_start_k(k);
        

        end

        
        % Check for the END of the current file- if it ends at the same time
        % step that another file begins OR if it is alone 

        % 2) Get x1 (input for GCC) - this sorts our the end of x array-
        % the last window of GCC

        % [L2]=ismember(indx_timestep_end_k(k),indx_timestep_start_k);
        if k<N
            [L2]=seconds(time_start(k+1) - global_time(indx_timestep_end_k(k)))<window_length_s;
        else %if k is the last file then there is no other file after, so L2 must be false
            L2=false;
        end

        if L2 % when there is another file starting at the same time step that k-th file ends
            file2_start_delay = seconds(time_start(k+1) - time_end(k)); % delay in start of the k+1 file relative to the k-th file (in seconds)


            if file2_start_delay>=0 % the files are non-overlapping
                %figure out how many samples are in the gap between x and x1 (so start of x1)

                missing_samples2=(file2_start_delay*fs); % delay in start of the file in samples
                 % missing_samples2=(file2_start_delay*fs)-1;

                %figure out how many samples to read to have a full window:
                end_time_step=find_upperdt(time_start(k+1),dt); %start of next time step in seconds
                end_sample= round((end_time_step-second(time_start(k+1))+window_length_s-dt)*fs);
                % (window_length_s-dt)*fs are missing sample to achieve a
                % full window for GCC.
                [x1,~]=audioread([folder1,'/',file1(k+1).name],[1,end_sample]);
                x1 = x1(:,channels);

                % zero pad the file
                x1 = [zeros(missing_samples2,size(x1,2));x1];

            else % the files are overlapping

                %figure out how many samples to read to have a full window:
                end_time_step=find_upperdt(time_start(k+1),dt);
                end_sample= round((end_time_step-second(time_start(k+1))+window_length_s-dt)*fs);
                
                start_sample= round(seconds((time_end(k)) - (time_start(k+1)))*fs);
                
                [x1,~]=audioread([folder1,'/',file1(k+1).name],[start_sample,end_sample]);
                x1 = x1(:,channels);

            end

        else % when there is no file starting at the same time step that k-th file ends
            x1=[];

        end


        stop_timestep=indx_timestep_end_k(k);
          
        %% Specify FILTERS and do FILTERING
        %~~~~~ Bandpass Butterworth ~~~~~~~~~
        [b,a] = butter(4,freq_filter/(fs/2),'bandpass'); %Coefficients for Butterworth
        x_filt = filtfilt(b,a,x); % Butterworth in forward-reverse fashion- so that it is zero phase
        x1_filt = filtfilt(b,a,x1);
        
        switch signal_type
            case 'whistles'
                %~~~~ Remove clicks and filter ~~~~~~~
                y = preprocess(x_filt,1024,0.5,5,6); %tested and it does not affect the phase
                y1 = preprocess(x1_filt,1024,0.5,5,6);
            case 'clicks'
                y=x_filt;
                y1=x1_filt;
            otherwise
                disp('Not a valid signal type. Choose "whistles" or "clicks"')
        end

        
        %% SPECIFY PARAMS for GCC and do GCC
        
        HT=1; %compute envelope via Hilbert transform
        
        [~,lags,RT_envelope,~] = gcc(y,y1,fs,window_length_samples,overlap,method,freq_filter,HT);

        ind = find(abs(lags)<parameters.tau_max); %get only lags that are physically possible
        lags=lags(ind);
        RT_envelope = RT_envelope(ind,:);

             
        Rxy_envelope_ALL(:,start_timestep:stop_timestep) = RT_envelope;


    end

    % %trim away unused pre-allocated space
    % Rxy_envelope_ALL=Rxy_envelope_ALL(1:m,:);

    t=(0:(size(Rxy_envelope_ALL,2)-1)).*dt; %start times of frames, in seconds (for easier plotting)

  
    %% SAVE WORKSPACE containing CROSS-CORELOGRAM
    
    if saveworksp==1
    switch signal_type
        case 'whistles'
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
        figure;
        imagesc(t_serialdate,lags,Rxy_envelope_ALL), datetick('x','keeplimits');
        colormap(flipud(gray(256)))
        colorbar
        ylim([-parameters.tau_max,parameters.tau_max])
        xlabel('Local Time (HH:MM)'), ylabel('TDOA (s)'),
        title(['Cross-correlogram based on ',signal_type])
        set(gca,'FontSize',12)
        set(findall(gcf,'type','text'),'FontSize',14)

    end

end

