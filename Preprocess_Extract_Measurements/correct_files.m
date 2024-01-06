function [x,Nstepsmissing] = correct_files(starttiming_mismatch,x,fs,dt)
% function that resolves missing or too many samples in the files. 
%
% INPUTS:
% - starttiming_mismatch - any potential mismatch between what was expected
%   in the start time of a file and what actually happened- this is in
%   seconds.
% - x - audiofile - Nx d matrix where N is number of samples and d is
%       number of channels.
% - fs - sample rate
% -dt - time step

if starttiming_mismatch<0 %If the file started too early, then chop off that data and discard
    chopoff_samples = abs(starttiming_mismatch)*fs;
    x = x(chopoff_samples+1:end,:);
    Nstepsmissing = 0;

elseif starttiming_mismatch>0 %If the file started too late either zero pad if it is less than one time
        %step or assign it to later (appropriate time steps).
        if starttiming_mismatch<dt
            missing_samples=starttiming_mismatch*fs;
            x = [zeros(missing_samples,size(x,2));x];
            Nstepsmissing =0;

        else %compute how many seconds are missing and convert to number of time steps
            missing_sec=starttiming_mismatch;
            Nstepsmissing=missing_sec/dt;

            %This is still problematic if I have for example a missing file
            %before this one, then this file starts prematurely- instead of
            %chopping off the premature part, the file stays as is- this
            %causes problems with GCC that then computes another frame..
            %RESOLVE
        end
else
    
    Nstepsmissing=0;
    
end

end