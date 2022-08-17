function [parameters,varargout]= Specify_Parameters(sigtype)
% SPECIFY PARAMETERS for ARRAY, ENCOUNTER &  CROSS-CORRELOGRAM:


if sum([strcmp(sigtype,'clicks'),strcmp(sigtype,'whistles'),strcmp(sigtype,'both')])==0
    msg= 'Not a valid signal type option. Choose "whistles" or "clicks" or "both".';
    error(msg)
end

%////////////////ARRAY and ENCOUNTER PARAMETERS & INFO/////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%Specify survey year
parameters.year=2017;

%Specify array name
parameters.arrayname= 'Numbat_20m_Goanna';

%Specify ENCOUNTER:
parameters.encounter= 'Lasker_AC276'; %sperm whale A196 uses numbat_gonna array so the same as Lasker AC191

% Specify speed of sound
parameters.c=1500; 

% Specify which two channels of your recordings you want to cross-correlate
% Select two that are furthest apart
parameters.channels = [1,6]; % Note All NOAA 2013 data has only four channels

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Compute distance between sensors used for cross-correlation:
Array_table=readtable('Array_Info.csv');
RowIn = find(strcmp(Array_table.Array_name, parameters.arrayname));
if ~isempty(RowIn)
    RowIn=RowIn(1);
    headers=Array_table.Properties.VariableNames;
    ColumnIn_1=contains(headers,num2str(parameters.channels(1)));
    ColumnIn_2=contains(headers,num2str(parameters.channels(2)));
    %Hydrophone separation:
    parameters.d=Array_table{RowIn,ColumnIn_2}-Array_table{RowIn,ColumnIn_1};
    if parameters.d<=0
        msg = ['Hydrophone distance error. Wrong hydrophone channels selected. ' ...
            ' Refer to Array_Info.csv to see which hydrophones are available. '];
        error(msg)
    end
else
    %     prompt="What is the sensor separation?";
    %     parameters.d=input(prompt);
    msg = ['The specified array is not found in Array_Info.csv ' ...
        '(folder Preprocess_Extract_Measureemnts). Please add the array ' ...
        'info into that file and try again.'];
    error(msg)
end
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////


%////////////////PARAMETERS for CROSS-CORRELOGRAM COMPUTATION/////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
% Specify the method for generalized cross-correlation (GCC):
%  'scc'/'phat'/'scot' (see gcc.m for more info)
parameters.method='scot'; %smoothed coherence transform

% Specify window length (in s) over which the cross-correlation is computed
parameters.window_length_s = 1; %1 s window

% Specify proportion of overlap between consecutive windows
parameters.overlap = 0.5; % 50% overlap

switch sigtype
    case 'clicks'
        varargout{1}.signal_type = 'clicks';
        % Bandpass Filter cut-off frequencies [lower, upper] (in Hz)
        % (determined based on which signal type is chosen):
        varargout{1}.freq_filter =[8000,30000]; % for false killer whales
        %         varargout{1}.freq_filter =[5500,20000]; %for sperm whales
%                 varargout{1}.freq_filter =[10000,80000]; %for Rissos dolphins
        parameters.signal_type = 'clicks'; %for plotting
    case 'whistles'
        varargout{1}.signal_type = 'whistles';
        varargout{1}.freq_filter =[2500,12000]; % for false killer whales
        parameters.signal_type = 'whistles'; %for plotting
    case 'both'
        varargout{1}.signal_type = 'clicks';
        % Bandpass Filter cut-off frequencies [lower, upper] (in Hz)
        % (determined based on which signal type is chosen):
        varargout{1}.freq_filter =[8000,30000]; % for false killer whales
        %         varargout{1}.freq_filter =[5500,20000]; %for sperm whales
        %         varargout{1}.freq_filter =[20000,80000]; %for Rissos dolphins
        varargout{2}.signal_type = 'whistles';
        varargout{2}.freq_filter =[2500,12000]; % for false killer whales
        %         varargout{2}.freq_filter =[2000,20000]; % for Rissos
        parameters.signal_type = 'both'; %for plotting
end

% Specify if you want the cross-correlogram figure plotted & saved (the
% cross-correlograms will also be plotted at the end in the plot_results.m
% function)
parameters.plotcrosscorr= 0; %Plot and save cross-correlogram (0=no,1=yes)

%Specify if you want the cross-correlogram saved separately. It will be
%saved as part of the RUN_TDOA_TRACKING script, so default is to not save
%separately.
parameters.saveworksp=0; % Save cross-correlogram & parameters (0=no,1=yes)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time increment between consecutive GCC windows (in s) (computed from 
% window length and overlap):
parameters.dt=(1-parameters.overlap)*parameters.window_length_s;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////


%////////////////PARAMETERS for MEASUREMENT EXTRACTION /////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%Specify threshold above which the cross-correlation peaks will be
%extracted:
parameters.lambda = sqrt(-2*log(0.001)); %Threshold lambda- we assume p_{FA}=0.001; Eq(5) Clark et al 2008

% Specify if you wish to exclude certain bearings in front of the array in
%order to avoid cross-correlation peaks due to the boat:
parameters.excl_deg = 10; % Cut off for bearings in degrees- e.g. 
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
parameters.excl_lags=-cosd(parameters.excl_deg)*parameters.d/parameters.c; 

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters.movingmeanlength=movingmeanlength_s/parameters.dt; %in time steps
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%//////////////////////////////////////////////////////////////////////////
end










