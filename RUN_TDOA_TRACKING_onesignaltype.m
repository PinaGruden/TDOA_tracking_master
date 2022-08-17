% RUN TDOA Tracking
% This is for one signal type- either clicks or whistles. If your species
% produces both clicks and whistles use RUN_TDOA_TRACKING_mixedsignaltypes.m.
% Also, make sure you specify and change parameters as needed in
% Specify_Paths.m and Specify_Parameters_CrossCorrelogram.m

clear,close all

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Choose 'whistles' or 'clicks
signal_type = 'clicks';

%Add current function and subfolders to the path:
% Determine where your m-file's folder is.
folder_code = fileparts(which('RUN_TDOA_TRACKING_onesignaltype.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder_code));

%% 1) USER SPECIFIED PATHS and PARAMETERS:

[folder_data,folder2save2] = Specify_Paths;

switch signal_type
    case 'clicks'%for clicks:
[parameters,parameters_clicks] = Specify_Parameters(signal_type);
    case 'whistles'% for whistles : 
[parameters,parameters_whistles] = Specify_Parameters(signal_type);
end

%% 2) CROSS-CORRELOGRAM COMPUTATION (and SAVE):

switch signal_type
    case 'clicks'%for clicks:
        [Rxy_envelope_ALL,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2, parameters,parameters_clicks);
    case 'whistles'% for whistles :
        [Rxy_envelope_ALL,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2, parameters,parameters_whistles);
end

%% 3) EXTRACTION OF MEASUREMENTS:

[measure,Rxy_envelope_scaled,scalar] = extract_measure_crosscorrelogram(...
    Rxy_envelope_ALL,lags,t, parameters.lambda, parameters.excl_lags,...
    parameters.tmax,parameters.min_dist);


%% 4) TRACK TDOAS:

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

plot_results(t_serialdate,lags,Rxy_envelope_ALL,measure,Tracks, parameters)

%% SAVE WORKSPACE with results &  measureemnts

switch signal_type
    case 'whistles'
        C=whos('Rxy_envelope_scaled');
        if C.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
            save([folder2save2,parameters.encounter,'_1s_05overlap_whistles_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled','scalar', 'lags','t_serialdate','t',...
                'parameters','measure','parameters_whistles', 'model', ...
                'Tracks','-v7.3')
        else
            save([folder2save2,parameters.encounter,'_1s_05overlap_whistles_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled','scalar', 'lags','t_serialdate','t',...
                'parameters','measure','parameters_whistles', 'model', ...
                'Tracks')
        end
    case 'clicks'
        C=whos('Rxy_envelope_scaled');
        if C.bytes>=2e+9 %if variable size is larger than 2GB then need to save as '-v7.3' (but it performs data compression so it's slow- so only use when necessary)
            save([folder2save2,parameters.encounter,'_1s_05overlap_clicks_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled','scalar', 'lags','t_serialdate','t',...
                'parameters','measure','parameters_clicks', 'model', ...
                'Tracks','-v7.3')
        else
            save([folder2save2,parameters.encounter,'_1s_05overlap_clicks_Rxy_envelope_ALL.mat'],...
                'Rxy_envelope_scaled','scalar', 'lags','t_serialdate','t',...
                'parameters','measure','parameters_clicks', 'model', ...
                'Tracks')
        end
    otherwise
        disp('Not a valid signal type. Choose "whistles" or "clicks"')
end
