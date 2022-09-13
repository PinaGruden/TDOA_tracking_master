%% COMPUTE CROSS-CORRELOGRAMS
% A1_Compute_CrossCorrelograms.m is a script that computes
% cross-correlograms based on two sensors. 
%!!!!!!!!!!!!!!!!!
% Make sure you first
% 1) Specify and change parameters as needed in
% Specify_Paths.m and Specify_Parameters4Xcorr.m
% 2) Specify signal_type variable to be either: 'clicks',
% 'whistles', or 'both'.
%!!!!!!!!!!!!!!!!!

addpath('./Preprocess_Extract_Measurements'); 
clear, close all


signal_type ='both';


%% 1) USER SPECIFIED PATHS and PARAMETERS:

[folder_data,folder2save2] = Specify_Paths;

switch signal_type
    case 'clicks'%for clicks:
        [parameters,parameters_clicks] = Specify_Parameters4Xcorr(signal_type);
    case 'whistles'% for whistles :
        [parameters,parameters_whistles] = Specify_Parameters4Xcorr(signal_type);
    case 'both'
        [parameters,parameters_clicks,parameters_whistles] = Specify_Parameters4Xcorr(signal_type);
end

parameters.saveworksp=1; %Save the overall cross-correlogram to folder2save2

%% 2) CROSS-CORRELOGRAM COMPUTATION (and SAVE)

disp('Computing cross-correlograms...')
switch signal_type
    case 'clicks'%for clicks:
        [Rxy_envelope_ALL,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2.rawcrosscorr, parameters,parameters_clicks);
    case 'whistles'% for whistles :
        [Rxy_envelope_ALL,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2.rawcrosscorr, parameters,parameters_whistles);
    case 'both'
        [Rxy_envelope_ALL_clicks]=compute_crosscorrelogram(folder_data,folder2save2.rawcrosscorr, parameters,parameters_clicks);
        [Rxy_envelope_ALL_whistles,lags,t_serialdate,t]=compute_crosscorrelogram(folder_data,folder2save2.rawcrosscorr, parameters,parameters_whistles);
    otherwise
        disp(['Not a valid signal type. signal_type needs to be either' ...
             ' "clicks","whistles", or "both". '])
end




