function [Lags,LagsAmp] = extract_peaks_crosscorr(Rxy,lags,lambda,varargin)
% extract_peaks_crosscorr.m is a function that extract peaks from the
% cross-correlogram. Optionally, one can specify the TDOAs to exclude (to
% avoid for example the vessel detection if towed array is used), and
% specify minimum distance between neighbouring peaks.
%
% INPUTS:
% - Rxy = Cross-correlogram- N x M matrix containing cross-correlation 
%   values, N number of TDOAs (lags), M number of time steps
% - lags = a 1 x N vector of TDOAs
% - lambda = a scalar specifying minimum amplitude threshold for peak selection 
% - varargin:
%   ~ varargin{1} = excl_lags = a value for TDOAs below which the search is
%       not performed (to avoid detecting a boat in front of the array)
%   ~ varargin{2} = min_dist = minimum distance between two neighbouring 
%       peaks of the cross-correlation function
%
% OUTPUTS:
% - Lags = extracted TDOAs- 1 x M cell array, where M is number of time steps 
% - LagsAmp = extracted amplitudes of the cross-correlation- 1 x M cell
%   array, where M is number of time steps
%
%
%Pina Gruden, UH Manoa, April 2022


switch nargin
    case 3 %default
        excl_lags = min(lags);% no lag gets excluded
        min_dist = 0.0031; % Deafult 0.0031- based on the simulation for
        %signals of 600 Hz bandwidth- if two sources are closer than that 
        %they will not be separated
    case 4
        excl_lags = varargin{1};
        min_dist = 0.0031; %Default
    case 5
        excl_lags = varargin{1};
        min_dist = varargin{2};
    otherwise
        error('This number of arguments is not supported')
end

%% Extract peaks
range=lags>excl_lags; %exclude the peaks due to boat

%EXTRACT PEAKS ABOVE THRESHOLD lambda:
Lags =cell(1,size(Rxy,2));
LagsAmp=Lags;

for k=1:size(Rxy,2) %Select all largest peaks above tau in each time step:
    warning('off', 'signal:findpeaks:largeMinPeakHeight') % to turn off warning of no peaks greater than MinPeakHeight
    [PkAmp, PkTime] = findpeaks(Rxy(range,k),lags(range),'MinPeakHeight',lambda,'MinPeakDistance',min_dist/2);
    Lags{k}=PkTime;
    LagsAmp{k}=PkAmp';
end



end