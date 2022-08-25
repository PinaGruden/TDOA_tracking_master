% RUN TDOA Tracking
% This is for one species that produce both signal types- clicks and whistles. 
% If your species produces only clicks or only whistles use 
% RUN_TDOA_TRACKING_onesignaltype.m.
% Also, make sure you specify and change parameters as needed in
% Specify_Paths.m and Specify_Parameters_CrossCorrelogram.m

clear,close all
signal_type = 'both';


%Add current function and subfolders to the path:
% Determine where your m-file's folder is.
folder_code = fileparts(which('RUN_TDOA_TRACKING_mixedsignaltypes.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder_code));

%% 1) USER SPECIFIED PATHS and PARAMETERS:

[folder_data,folder2save2] = Specify_Paths;

[parameters,parameters_clicks,parameters_whistles] = Specify_Parameters(signal_type);

%% 2) CROSS-CORRELOGRAM COMPUTATION (and SAVE):

[Rxy_envelope_ALL_clicks]=compute_crosscorrelogram(folder_data,folder2save2, parameters,parameters_clicks);
[Rxy_envelope_ALL_whistles,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2, parameters,parameters_whistles);


%% 3) EXTRACTION OF MEASUREMENTS:
disp('Extracting measurements')
[measure_clicks,Rxy_envelope_scaled_clicks,scalar_clicks] = extract_measure_crosscorrelogram(...
    Rxy_envelope_ALL_clicks,lags,t, parameters.lambda, parameters.excl_lags,...
    parameters.tmax,parameters.min_dist);
[measure_whistles,Rxy_envelope_scaled_whistles,scalar_whistles] = extract_measure_crosscorrelogram(...
    Rxy_envelope_ALL_whistles,lags,t, parameters.lambda, parameters.excl_lags,...
    parameters.tmax,parameters.min_dist);

%merge measurements:
measure.T=measure_clicks.T;
measure.Z= cellfun(@(x,y) [x,y],measure_clicks.Z,measure_whistles.Z,'UniformOutput',false);
%% 4) TRACK TDOAS:
disp('Tracking TDOAs from measurements')
% Generate tracking models:
load('BayesOptimization_GMPHDParams_GTchunked_lambda4_5_DiffBirth.mat')
load('BirthVelocityPrior_AllTrainData.mat')
[model] = gen_tdoa_tracking_models(parameters, BayesoptResults,birthvelocity);

% Employ GM-PHD-SA to track:
Est = gmphd_adaptive_amplitude(model,measure); 
Track = tracktarget_tdoa_labels(Est.Tag,Est.X,model);

% Interpolate and smooth tracks, apply min track length criteria:
[Tracks] = postprocess(Track,t_serialdate,parameters);


%% 5) PLOT:
disp('Plotting results')
plot_results(t_serialdate,lags,{Rxy_envelope_scaled_clicks,Rxy_envelope_scaled_whistles},measure,Tracks, parameters)

%% 6) SAVE WORKSPACE with results &  measureemnts

        Cc=whos('Rxy_envelope_scaled_clicks');
        Cw=whos('Rxy_envelope_scaled_whistles');
        if Cc.bytes>=2e+9 || Cw.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
            save([folder2save2,parameters.encounter,'_1s_05overlap_both_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled_clicks','scalar_clicks', ...
                'Rxy_envelope_scaled_whistles','scalar_whistles', ...
                'lags','t_serialdate','t',...
                'parameters','parameters_whistles', 'parameters_clicks', ...
                'measure','model','Tracks','-v7.3')
        else
            save([folder2save2,parameters.encounter,'_1s_05overlap_both_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled_clicks','scalar_clicks', ...
                'Rxy_envelope_scaled_whistles','scalar_whistles', ...
                'lags','t_serialdate','t',...
                'parameters','parameters_whistles', 'parameters_clicks', ...
                'measure','model','Tracks')
        end
        
    