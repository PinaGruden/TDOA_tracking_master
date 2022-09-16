function [parameters]= Specify_Parameters4Tracking(sigtype,paramsXcorr)
% SPECIFY PARAMETERS for MEASUREMENT EXTRACTION and TRACKING:


if sum([strcmp(sigtype,'clicks'),strcmp(sigtype,'whistles'),strcmp(sigtype,'both')])==0
    msg= 'Not a valid signal type option. Choose "whistles" or "clicks" or "both".';
    error(msg)
end


%////////////////PARAMETERS for MEASUREMENT EXTRACTION /////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%Specify threshold above which the cross-correlation peaks will be
%extracted:
parameters.lambda = sqrt(-2*log(0.001)); %Threshold lambda- we assume p_{FA}=0.001; Eq(5) Clark et al 2008

% Specify if you wish to exclude certain bearings in front of the array in
%order to avoid cross-correlation peaks due to the boat:
parameters.excl_deg = 15; % Cut off for bearings in degrees- e.g. 
% no bering smaller than 10 degrees will be considered

% Specify maximum time (in s) to be considered in the noise sample- 
%   time from the beggining of the encounter up to before any significant 
%   sources occur (e.g. 1000s).
parameters.tmax=100;

% Specify minimum distance between two neighbouring cross-correlation peaks
parameters.min_dist=0.0031; % For FKW whistles- based on the simulation for
        %signals of 600 Hz bandwidth- if two sources are closer than that 
        %they will not be separated

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%A threhold below which TDOAs are not considered - exclude TDOAs directly in front of the array- boat 
%( the way we specify the TDOAs ahead are negative)
parameters.excl_lags=-cosd(parameters.excl_deg)*paramsXcorr.d/paramsXcorr.c; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%//////////////////////////////////////////////////////////////////////////

%////////////////PARAMETERS for EXTRACTED TRACKS ///////////////////////
% This is used in post-processing to interpolate and smooth GMPHD-SA tracks

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
% Specify minimum track length criteria- a TDOA track 
% must be at least min_tl long in order to become a true track:
parameters.min_tl=3; %in s

% Specify the length of the moving average (for track smoothing) - if no
% track smoothing & interpolation is desired then set movingmeanlength_s=0;
movingmeanlength_s=5; %in s

%Save result Plots - specify if you want to save plots of your results (1=yes,0 = no)
parameters.saveplotsofresults=1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters.movingmeanlength=movingmeanlength_s/paramsXcorr.dt; %in time steps
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%//////////////////////////////////////////////////////////////////////////
end









