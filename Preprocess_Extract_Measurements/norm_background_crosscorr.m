function [Rxy_scaled,scalar] = norm_background_crosscorr(Rxy,lags,t,varargin)
%norm_background is a function to normalize background noise to have std =
%1 under the Rayleigh assumption. The assumption is that the negative lags
%refer to positions ahead of the array, and positive lags are behind the
%array.
%
% INPUTS:
% -Rxy= Cross-correlogram - N x M matrix containing cross-correlation 
%         values, N= length(lags), M = length(t)
% - lags = a 1 x N vector of TDOAs,
% - t = a 1 x M vector of times (in s)
% -varargin:
%   ~ varargin{1} = excl_lags= a threshold value for TDOAs below which the 
%       TDOAs are excluded from further analysis (to avoid detecting a boat  
%       in front of the array)
%   ~ varargin{2} = tmax = maximum time (s) to be considered in the noise  
%       sample- time from the beggining of the encounter up to before any 
%       significant sources occur
%
% OUTPUTS:
% - Rxy_scaled= Rayleigh normalized cross-correlogram (Rayleigh parameter 
%   sigma = 1) - N x M matrix containing normalized cross-correlation values,
%    where N is number of TDOAs, M is number of time steps.
% - scalar = scalar used to normalize the cross-correlogram
%
%
%Pina Gruden, UH Manoa, 2022

switch nargin
    case 3 %default
        excl_lags= min(lags); %default includes all TDOA values
        tmax=max(t); %default includes all time steps
    case 4
       excl_lags= varargin{1};
       tmax=max(t); %default includes all time steps
    case 5
        excl_lags= varargin{1};
        tmax= varargin{2};
    otherwise
        error('This number of arguments is not supported')
end

 
Rxy_envelope_noiseonly=Rxy(lags>excl_lags,t<=tmax);
Rxy_envelope_noiseonly=Rxy_envelope_noiseonly(:);

[sigma,~] = raylfit(Rxy_envelope_noiseonly); %obtain the std
scalar = 1/sigma; %this ensure that the Rayleigh parameter for noise case will be 1 (sigma=1)
Rxy_scaled = Rxy.*scalar;

end