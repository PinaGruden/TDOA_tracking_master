%% EXTRACT MEASUREMENTS
% A2_Extract_Measurements.m is a script that extracts measurements from
% cross-correlograms. 
%!!!!!!!!!!!!!!!!!
% IMPORTANT: This assumes that data for one encounter at the time is being
% processed.
% Make sure you have first
% 1) Already run A1_Compute_CrossCorrelograms.m & specified paths in
% Specify_paths.m (should have been specified before running A1)
%!!!!!!!!!!!!!!!!!

addpath('./Preprocess_Extract_Measurements'); 
clear, close all

[~,folder2save2] = Specify_Paths;

folder_crosscorr = folder2save2.rawcrosscorr;


file = dir(fullfile(folder_crosscorr,'*.mat' ) ); % list all wav files in subfolder
%-----------------------------------------------------------------
%eliminate all '.' and '..' files:
%-----------------------------------------------------------------
indx=[];
for k = 1:size(file,1) %iterate through files in the folder
    if file(k).name(1) == '.'
        indx=[indx,k]; % skip '.' and '..' files
    end
end
file(indx)=[];
clearvars indx k 
N=size(file,1);


%% EXTRACT MEASUREMENTS (and SAVE):

if N==1 % Only Clicks or Whistles were used to create Cross-correlograms
    % LOAD CROSS-CORRELOGRAM
    load([file.folder,'/',file.name]);

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS
    disp('Extracting measurements from cross-correlogram...')
    [measure,Rxy_envelope_scaled,scalar] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters.lambda, parameters.excl_lags,...
        parameters.tmax,parameters.min_dist);

    % SAVE
    C=whos('Rxy_envelope_scaled');
    if C.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
        save([folder2save2.measurements,parameters.encounter,'_Measurements.mat'],...
            'measure','Rxy_envelope_scaled', 'scalar', ...
            'lags','t_serialdate', 't',...
            'parameters','param_signal','-v7.3')
    else
        save([folder2save2.measurements,parameters.encounter,'_Measurements.mat'],...
            'measure','Rxy_envelope_scaled', 'scalar', ...
            'lags','t_serialdate', 't',...
            'parameters','param_signal')
    end

elseif N==2 %Both Clicks and Whistles were used to create Cross-correlograms
    % 1) CLICKS:
    % LOAD CROSS-CORRELOGRAM
    load([file(1).folder,'/',file(1).name]);
    parameters_clicks=param_signal;

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS:
    disp('Extracting measurements from click cross-correlogram...')
    [measure_clicks,Rxy_envelope_scaled_clicks,scalar_clicks] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters.lambda, parameters.excl_lags,...
        parameters.tmax,parameters.min_dist);

    % 2) WHISTLES:
    % LOAD CROSS-CORRELOGRAM
    load([file(2).folder,'/',file(2).name]);
    parameters_whistles=param_signal;

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS:
    disp('Extracting measurements from whistle cross-correlogram...')
    [measure_whistles,Rxy_envelope_scaled_whistles,scalar_whistles] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters.lambda, parameters.excl_lags,...
        parameters.tmax,parameters.min_dist);

    % 3) MERGE measurements:
    measure.T=measure_clicks.T;
    measure.Z= cellfun(@(x,y) [x,y],measure_clicks.Z,measure_whistles.Z,'UniformOutput',false);


    % SAVE
    Cc=whos('Rxy_envelope_scaled_clicks');
    Cw=whos('Rxy_envelope_scaled_whistles');
    if Cc.bytes>=2e+9 || Cw.bytes>=2e+9%if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
        save([folder2save2.measurements,parameters.encounter,'_Measurements.mat'],...
            'measure','Rxy_envelope_scaled_clicks','scalar_clicks', ...
            'Rxy_envelope_scaled_whistles','scalar_whistles', ...
            'lags','t_serialdate', 't',...
            'parameters','parameters_whistles', 'parameters_clicks','-v7.3')
    else
        save([folder2save2.measurements,parameters.encounter,'_Measurements.mat'],...
            'measure','Rxy_envelope_scaled_clicks','scalar_clicks', ...
            'Rxy_envelope_scaled_whistles','scalar_whistles', ...
            'lags','t_serialdate', 't',...
            'parameters','parameters_whistles','parameters_clicks')
    end


else
    disp(['Wrong number of cross-correlograms. This processing scheme assumes ' ...
        '1 cross-correlogram per encounter for clicks and/or ' ...
        '1 cross-correlogram per encounter for whistles. Make ' ...
        'sure your raw cross-correlogram folder only contains data from ' ...
        'one encounter'])
end






