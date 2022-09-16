%% EXTRACT MEASUREMENTS and TRACK TDOAs
% A2_ExtractMeasurements_and_TrackTDOAs.m is a script that extracts measurements from
% cross-correlograms. 
%!!!!!!!!!!!!!!!!!
% IMPORTANT: This assumes that data for one encounter at the time is being
% processed.
% Make sure you have first
% 1) Already run A1_Compute_CrossCorrelograms.m & specified paths in
% Specify_paths.m (should have been specified before running A1)
% 2) Specified parameters that you want to use for measurement extraction
% and tracking in Specify_Parameters4Tracking.m
%!!!!!!!!!!!!!!!!!

addpath('./Preprocess_Extract_Measurements'); 
addpath('./GMPHD_SA'); 
clear, close all

[~,folder2save2] = Specify_Paths;

folder_crosscorr = folder2save2.rawcrosscorr;


file = dir(fullfile(folder_crosscorr,'*.mat' ) ); % list all .mat files in subfolder
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


%% EXTRACT MEASUREMENTS:


if N==1 % Only Clicks or Whistles were used to create Cross-correlograms
    % LOAD CROSS-CORRELOGRAM
    load([file.folder,'/',file.name]);

    % Get parameters for measurements & tracking:
    [parameters_measure_tracking]= Specify_Parameters4Tracking(parameters.signal_type,parameters);

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS
    disp('Extracting measurements from cross-correlogram...')
    [measure,Rxy_envelope_scaled,scalar] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters_measure_tracking.lambda, parameters_measure_tracking.excl_lags,...
        parameters_measure_tracking.tmax,parameters_measure_tracking.min_dist);

    
elseif N==2 %Both Clicks and Whistles were used to create Cross-correlograms
    % 1) CLICKS:
    % LOAD CROSS-CORRELOGRAM
    load([file(1).folder,'/',file(1).name]);
    parameters_clicks=param_signal;

    % Get parameters for measurements & tracking (this wil be the same for clicks and whistles):
    [parameters_measure_tracking]= Specify_Parameters4Tracking(parameters.signal_type,parameters);

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS:
    disp('Extracting measurements from click cross-correlogram...')
    [measure_clicks,Rxy_envelope_scaled_clicks,scalar_clicks] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters_measure_tracking.lambda, parameters_measure_tracking.excl_lags,...
       parameters_measure_tracking.tmax,parameters_measure_tracking.min_dist);

    % 2) WHISTLES:
    % LOAD CROSS-CORRELOGRAM
    load([file(2).folder,'/',file(2).name]);
    parameters_whistles=param_signal;

    % NORMALIZE CROSS-CORRELOGRAM & EXTRACT MEASUREMENTS:
    disp('Extracting measurements from whistle cross-correlogram...')
    [measure_whistles,Rxy_envelope_scaled_whistles,scalar_whistles] = extract_measure_crosscorrelogram(...
        Rxy_envelope_ALL,lags,t, parameters_measure_tracking.lambda, parameters_measure_tracking.excl_lags,...
        parameters_measure_tracking.tmax,parameters_measure_tracking.min_dist);

    % 3) MERGE measurements:
    measure.T=measure_clicks.T;
    measure.Z= cellfun(@(x,y) [x,y],measure_clicks.Z,measure_whistles.Z,'UniformOutput',false);


else
    disp(['Wrong number of cross-correlograms. This processing scheme assumes ' ...
        '1 cross-correlogram per encounter for clicks and/or ' ...
        '1 cross-correlogram per encounter for whistles. Make ' ...
        'sure your raw cross-correlogram folder only contains data from ' ...
        'one encounter'])
end


%% TRACK TDOAs:

parameters.lambda=parameters_measure_tracking.lambda;
parameters.min_tl=parameters_measure_tracking.min_tl;
parameters.movingmeanlength=parameters_measure_tracking.movingmeanlength;

disp('Tracking TDOAs from measurements...')
% Generate tracking models:
load('BayesOptimization_GMPHDParams_GTchunked_lambda4_5_DiffBirth.mat')
load('BirthVelocityPrior_AllTrainData.mat')

[model] = gen_tdoa_tracking_models(parameters, BayesoptResults,birthvelocity);

% Employ GM-PHD-SA to track:
Est = gmphd_adaptive_amplitude(model,measure); 
Track = tracktarget_tdoa_labels(Est.Tag,Est.X,model);

% Interpolate and smooth tracks, apply min track length criteria:
[Tracks] = postprocess(Track,t_serialdate,parameters);


%% PLOT RESULTS:
disp('Plotting results...')
parameters.saveplotsofresults=parameters_measure_tracking.saveplotsofresults;

switch parameters.signal_type
    case 'clicks' 
        plot_results(t_serialdate,lags,Rxy_envelope_scaled,measure,Tracks, parameters,folder2save2.finalresults);
    case 'whistles'
        plot_results(t_serialdate,lags,Rxy_envelope_scaled,measure,Tracks, parameters,folder2save2.finalresults);
    case 'both'
        plot_results(t_serialdate,lags,{Rxy_envelope_scaled_clicks,Rxy_envelope_scaled_whistles},measure,Tracks, parameters,folder2save2.finalresults);
end

%% SAVE:

if N==1 % Only Clicks or Whistles were used to create Cross-correlograms

    save([folder2save2.finalresults,parameters.encounter,'_Results.mat'],...
        'model','Tracks','scalar','parameters','param_signal','parameters_measure_tracking')

else  %Both Clicks and Whistles were used to create Cross-correlograms

    save([folder2save2.finalresults,parameters.encounter,'_Results.mat'],...
        'model','Tracks','scalar_clicks', 'scalar_whistles', ...
            'parameters','parameters_whistles','parameters_clicks','parameters_measure_tracking')

end