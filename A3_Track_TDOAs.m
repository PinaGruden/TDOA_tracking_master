%% TRACK TDOAs and DISPLAY RESULTS
% A3_Track_TDOAs.m is a script that tracks TDOAs based on extracted 
% measurements from cross-correlograms. 
%!!!!!!!!!!!!!!!!!
% IMPORTANT: This assumes that data for one encounter at the time is being
% processed. This means only one .mat file that contains measurements for
% the full encounter is in a given folder.
% Make sure you have first
% 1) Already run A1_Compute_CrossCorrelograms.m & specified paths in
% Specify_paths.m (should have been specified before running A1)
% 2) Already run A2_Extract_Measurements.m
%!!!!!!!!!!!!!!!!!

addpath('./GMPHD_SA'); 
clear, close all

[~,folder2save2] = Specify_Paths;

file = dir(fullfile(folder2save2.measurements,'*.mat' ) );
load([file.folder,'/',file.name]);

%% TRACK TDOAs:
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

switch parameters.signal_type
    case 'clicks' 
        plot_results(t_serialdate,lags,Rxy_envelope_scaled,measure,Tracks, parameters);
    case 'whistles'
        plot_results(t_serialdate,lags,Rxy_envelope_scaled,measure,Tracks, parameters);
    case 'both'
        plot_results(t_serialdate,lags,{Rxy_envelope_scaled_clicks,Rxy_envelope_scaled_whistles},measure,Tracks, parameters);
end

%% SAVE

save([folder2save2.finalresults,parameters.encounter,'_Results.mat'],...
            'model','Tracks')