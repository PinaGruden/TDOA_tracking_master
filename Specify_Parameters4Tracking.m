function [parameters]= Specify_Parameters4Tracking(sigtype,paramsXcorr)
% Use Specify_Parameters4Tracking.m to SPECIFY PARAMETERS for MEASUREMENT 
% EXTRACTION and TRACKING. Change all parameters in sections labeled
% "CHANGABLE".
%
% INPUTS:
% - sigtype - a string specifying a signal type that is of interest. Needs
%             to be 'whistles' or 'clicks' or 'both'.
% - paramsXcorr - a structure with at least 3 fields:
%       ~ d - a number specifying sensor separation (m)
%       ~ c - a number specifying the speed of sound (m/s)
%       ~ dt - a number specifying time increment between consecutive time 
%       steps (in s)
%
% OUTPUTS:
% - parameters - a structure with 8 fields:
%       ~ lambda - a number specifying threshold above which the 
%           cross-correlation peaks will be extracted
%       ~ excl_deg - a number specifying a cut-off for bearings in degrees
%       ~ tmax - maximum time (in s) to be considered in the noise sample
%       ~ min_dist - a number specifying minimum distance between two 
%           neighbouring peaks in the  cross-correlation function
%       ~ excl_lags - a number specifying a threhold below which TDOAs are 
%           not considered
%       ~ min_tl - a number specifying a minimum track length criteria 
%           (in s)
%       ~ saveplotsofresults - a number specifying whether to save plots of
%           results (1=yes,0 = no)
%       ~ movingmeanlength - a number specifying the length of the moving 
%           average (for track smoothing) (in time steps)
%
%
% Pina Gruden, 2022, UH Manoa

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
% time from the beginning of the encounter until just before any 
% significant sources occur (e.g. 1000s). Note, the longer the sample the
% better the background estimation.
parameters.tmax=100;

% Specify minimum distance (in seconds) between two neighbouring cross-correlation peaks
%!!!!!!!!!!!!!!!
% NOTE: minumum distance (in seconds) should be LESS than 2 * max possible
%TDOA (determined by your sensor separation- 2*max_{tdoa} = 2*d/c -
%specified in Specify_Parameters4Xcorr.m)
%!!!!!!!!!!!!!!!!
parameters.min_dist=0.0031; % For FKW whistles- based on the simulation for
        %signals of 600 Hz bandwidth- if two sources are closer than that 
        %they will not be separated

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Check if parameters.min_dist < paramsXcorr.tau_max and if not warn the user
% to either decrease min_dist OR increase sensor separation d in
% Specify_Parameters4Xcorr.m - i.e. use different channels.
if parameters.min_dist>=paramsXcorr.tau_max
    error(['The chosen parameters.min_dist- i.e. minimum distance between two ' ...
        'neighboring peaks in cross-correlation is BIGGER than what is ' ...
        'physically possible - i.e. it exceeds 2 x maximum possible TDOA ' ...
        'determined based on sensor separation d in Specify_Parameters4Xcorr.m.' ...
        ' EITHER decrease the min_dist (in Specify_Parameters4Tracking.m) ' ...
        'OR choose to process sensors with bigger separation (in Specify_Parameters4Xcorr.m).'])
end

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










