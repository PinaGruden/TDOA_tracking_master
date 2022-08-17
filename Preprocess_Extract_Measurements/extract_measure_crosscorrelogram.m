function [measure,Rxy_envelope_scaled,scalar] = extract_measure_crosscorrelogram(Rxy_envelope_ALL,lags,t, lambda, varargin)
%extract_measure_crosscorrelogram.m is a function that extracts TDOA and
%amplitude of the cross-correlation information form a cross-correlogram
%and saves it as a random finite set (RFS).
%
% INPUTS:
% - RT_envelope_ALL = Cross-correlogram- N x M matrix containing 
%   cross-correlation values, N= length(lags), M = length(t)
% - lags = TDOAs,
% - t = time vector (in s)
%varargin:
% - varargin{1}= excl_lags = a value for TDOAs below which the TDOAs are
%   excluded from further analysis- this is to avoid detecting a boat in 
%   front of the array (negative TDOAs). Used in
%   norm_background_crosscorr.m and extract_peaks_crosscorr.m
% - varargin{2}= tmax = maximum time (s) to be considered in the noise sample- 
%   time from the beggining of the encounter up to before any significant 
%   sources occur (e.g. 1000s). Used in norm_background_crosscorr.m
% - varargin{3}= min_dist = minimum distance between two neighbouring
%   peaks. Used in extract_peaks_crosscorr.m 
%
% OUTPUTS:
%- measure = RFS of measurements containing TDOA and amplitude 
%   of the cross-correlation information
%- RT_envelope_scaled = Rayleigh normalized cross-correlogram,
%- scalar = scalar used to normalize cross-correlogram


%Pina Gruden, UH Manoa, April 2022


switch nargin
    case 4
        % \\\\\\\\\\\\\\\\\\\\ NORMALIZE AMPLITUDE DATA\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [Rxy_envelope_scaled,scalar] = norm_background_crosscorr(Rxy_envelope_ALL,lags,t);
        % ///////////////// EXTRACT TDOA PEAKS & AMPLITUDES /////////////////////
        [Lags,LagsAmp] = extract_peaks_crosscorr(Rxy_envelope_scaled,lags,lambda);
    case 5
        % \\\\\\\\\\\\\\\\\\\\ NORMALIZE AMPLITUDE DATA\\\\\\\\\\\\\\\\\\\\\\\\\\\
        excl_lags = varargin{1};
        [Rxy_envelope_scaled,scalar] = norm_background_crosscorr(Rxy_envelope_ALL,lags,t,excl_lags);
        % ///////////////// EXTRACT TDOA PEAKS & AMPLITUDES /////////////////////
        [Lags,LagsAmp] = extract_peaks_crosscorr(Rxy_envelope_scaled,lags,lambda,excl_lags);
    case 6
        % \\\\\\\\\\\\\\\\\\\\ NORMALIZE AMPLITUDE DATA\\\\\\\\\\\\\\\\\\\\\\\\\\\
        excl_lags = varargin{1};
        tmax=varargin{2};
        [Rxy_envelope_scaled,scalar] = norm_background_crosscorr(Rxy_envelope_ALL,lags,t,excl_lags,tmax);
        % ///////////////// EXTRACT TDOA PEAKS & AMPLITUDES /////////////////////
        [Lags,LagsAmp] = extract_peaks_crosscorr(Rxy_envelope_scaled,lags,lambda,excl_lags);
    case 7
        % \\\\\\\\\\\\\\\\\\\\ NORMALIZE AMPLITUDE DATA\\\\\\\\\\\\\\\\\\\\\\\\\\\
        excl_lags = varargin{1};
        tmax=varargin{2};
        [Rxy_envelope_scaled,scalar] = norm_background_crosscorr(Rxy_envelope_ALL,lags,t,excl_lags,tmax);
        % ///////////////// EXTRACT TDOA PEAKS & AMPLITUDES //////////////////////
        min_dist=varargin{3};
        [Lags,LagsAmp] = extract_peaks_crosscorr(Rxy_envelope_scaled,lags,lambda,excl_lags,min_dist);
    otherwise
        error('This number of arguments is not supported')
end


%% //////////////////// CREATE RFS MEASUREMENT SET //////////////////////////////
%Combine TDOAS and AMPLITUDES from the same time step:
N=size(Rxy_envelope_scaled,2);
measure.Z = cell(1,N);
for k=1:N
    measure.Z{k}=[Lags{k};LagsAmp{k}];
end
measure.T=N;

end