function [parameters,varargout]= Specify_Parameters4Xcorr(sigtype)
% SPECIFY PARAMETERS for ARRAY, ENCOUNTER &  CROSS-CORRELOGRAM:


if sum([strcmp(sigtype,'clicks'),strcmp(sigtype,'whistles'),strcmp(sigtype,'both')])==0
    msg= 'Not a valid signal type option. Choose "whistles" or "clicks" or "both".';
    error(msg)
end

%////////////////ARRAY and ENCOUNTER PARAMETERS & INFO/////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%Specify survey year
parameters.year=2022;

%Specify array name
parameters.arrayname= 'test';

%Specify ENCOUNTER:
parameters.encounter= 'Test'; 

% Specify speed of sound
parameters.c=1500; 

% Specify which two channels of your recordings you want to cross-correlate
% Select two that are furthest apart
parameters.channels = [1,2]; % Note All NOAA 2013 data has only four channels

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
        '(found in the main folder). Please add the array ' ...
        'info to it and try again.'];
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
% varargout{1}.freq_filter =[10000,50000]; % for rough toothed dolphins (Rankin et al 2008)
        parameters.signal_type = 'clicks'; %for plotting
    case 'whistles'
        varargout{1}.signal_type = 'whistles';
        varargout{1}.freq_filter =[2500,12000]; % for false killer whales
%         varargout{1}.freq_filter =[3300,28000]; % for rough toothed dolphins (Rankin et al 2008)
        parameters.signal_type = 'whistles'; %for plotting
    case 'both'
        varargout{1}.signal_type = 'clicks';
        % Bandpass Filter cut-off frequencies [lower, upper] (in Hz)
        % (determined based on which signal type is chosen):
        varargout{1}.freq_filter =[8000,30000]; % for false killer whales
        %         varargout{1}.freq_filter =[5500,20000]; %for sperm whales
        %         varargout{1}.freq_filter =[20000,80000]; %for Rissos dolphins
%         varargout{1}.freq_filter =[10000,50000]; % for rough toothed dolphins (Rankin et al 2008)
        varargout{2}.signal_type = 'whistles';
        varargout{2}.freq_filter =[2500,12000]; % for false killer whales
        %         varargout{2}.freq_filter =[2000,20000]; % for Rissos
%         varargout{2}.freq_filter =[3300,28000]; % for rough toothed dolphins (Rankin et al 2008)
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



end










