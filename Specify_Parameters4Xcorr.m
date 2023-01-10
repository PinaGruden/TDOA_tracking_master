function [parameters,varargout]= Specify_Parameters4Xcorr(sigtype)
% Use Specify_Parameters4Xcorr.m to SPECIFY PARAMETERS for ARRAY, ENCOUNTER 
% & CROSS-CORRELOGRAM. Change all parameters in sections labeled
% "CHANGABLE".

% INPUTS:
% - sigtype - a string specifying a signal type that is of interest. Needs
%             to be 'whistles' or 'clicks' or 'both'.

% OUTPUTS:
% - parameters - a structure containing parameters for the encounter and 
%               processing. It has 13 fields:
%       ~ parameters.year - a scalar specifying the year of the survey
%       ~ parameters.arrayname - a string specifying the array name (should
%         match the name in Array_info.csv)
%       ~ parameters.encounter- a string specifying the encounter name (should
%         match the name in Array_info.csv)
%       ~ parameters.c - a scalar specifying the speed of sound
%       ~ parameters.channels - 1 x 2 vector specifying sensor channels to
%           be used for processing
%       ~ parameters.d - a scalar specifying sensor separation
%       ~ parameters.method - a string specifying method for generalized
%           cross-correlation (GCC) computation. Should be 'scc'/'phat'/'scot' 
%           (see gcc.m for more info)
%       ~ parameters.window_length_s - window length (in s) over which the 
%           cross-correlation is computed
%       ~ parameters.overlap - proportion of overlap between consecutive
%           windows for cross-correlogram computation (between 0 and 1)
%       ~ parameters.signal_type - a string specifying a signal type that 
%           is of interest ('whistles'/'clicks'/'both')
%       ~ parameters.plotcrosscorr- a scalar to specify to plot the 
%           cross-correlogram or not (0=no,1=yes)
%       ~ parameters.saveworksp - a scalar to specify to save the 
%           cross-correlogram or not (0=no,1=yes)
%       ~ parameters.dt - time increment between consecutive GCC windows (in s)
% - vargout - a structure with two fields:
%       ~ signal_type -  a string specifying the signal type ('whistles'/'clicks') 
%       ~ freq_filter - 1 x 2 vector specifying min and max frequency limits
%        for the signal_type (in Hz).


% Pina Gruden, 2022, UH Manoa


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
parameters.channels = [1,2]; 

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

% Specify if you want the cross-correlogram figure plotted after the 
% cross-correlogram is computed - this can help troubleshoot (the
% cross-correlograms will also be plotted at the end after tracking 
% in A2_ExtractMeasurements_and_TrackTDOAs.m)
parameters.plotcrosscorr= 1; %Plot cross-correlogram (0=no,1=yes)

%Specify if you want the cross-correlogram saved.
parameters.saveworksp=1; % Save cross-correlogram & parameters (0=no,1=yes)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time increment between consecutive GCC windows (in s) (computed from 
% window length and overlap):
parameters.dt=(1-parameters.overlap)*parameters.window_length_s;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////



end










