function [Rxy_scaled,scalar] = norm_background_crosscorr(Rxy,lags,t,varargin)
%norm_background is a function to normalize background noise to have std =
%1 under the Rayleigh assumption. The assumption is that the negative lags
%refer to positions ahead of the array, and positive lags are behind the
%array.
%
%   Inputs:
%   -Rxy= cross-correlogram- m x n dimensional matrix,
%   -lags= 1 x m vector of TDOAs (lags) corresponding to the cross-correlogram,
%   -t= 1 x n vector of times (in s) corresponding to the cross-correlogram,
%   -varargin:
%       ~ varargin{1} = excl_lags= a threhold below which TDOAs are not considered -
%        e.g. if there is boat source ahead of the array exclude these values,
%       ~ varargin{2} = tmax = maximum time to be considered- time from the beggining of the
%   encounter up to before any significant sources occur.
%
%   Outputs:
%   -Rxy_scaled= scaled cross-correlogram to have a normalized background
%   noise with Rayleigh parameter sigma = 1.
%   -scalar = scalar used to normalize the cross-correlogram


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
Rxy_scaled = Rxy*scalar;

end