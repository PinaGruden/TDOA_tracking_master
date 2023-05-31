function [tracks] = postprocess(Track,t_serialdate,parameters)
% function postprocess.m removes short tracks (below min length criteria)
% form a list of tracks, smoothes the tracks trajectories (with moving
% average) and adds time info (in a serial date format) to each track.
%
% INPUTS:
% - Track - a structure contatining each target (TDOA track) in a separate
%       row. It has 5 fields:
%       ~ tdoa - tdoas of a given target track
%       ~ dottdoa - derivatives of the tdoas (i.e. velocity) of a given target track
%       ~ time - times (from the begining of the encounter) of a given target track
%       ~ label - labels (identities) of a given target track
%       ~ ti - index to a time step (from the beginning of the encounter) of a given target track
% - t_serialdate - a 1 x M vector of times (in serial date format) for the
%       encounter
% - parameters - a structure containing info on encounter, array, and
%       parameters used for measurement extraction. Needs at least 3 fields:
%       ~ dt - time step (in s)
%       ~ min_tl - minimum track length criteria (in s)
%       ~ movemean_length- length of the moving average (in time steps)
%
% OUTPUTS:
% - tracks - a structure contatining each target (TDOA track) in a separate
%       row. It has 3 fields:
%       ~ time - times (from the begining of the encounter) of a given target track
%       ~ time_local - times (in a serial date format) of a given target track
%       ~ tdoa - tdoas of a given target track
%
%
% Pina Gruden

dt=parameters.dt; %time step 
min_tl= parameters.min_tl; %minimum track length criteria (in s)
movemean_length = parameters.movingmeanlength; %length of the moving average


%Remove all tracks that have only one estimate:
N=size(Track,2);
indx=false(N,1);
for k=1:N
    if numel(Track(k).time)>1
        indx(k)=true;
    end
end
Track=Track(indx);

%interpolate to remove holes in individual track, and smooth
N=size(Track,2);
tracks(N).time=[];
tracks(N).time_local=[];
tracks(N).ti=[];
tracks(N).tdoa=[];

for k=1:N

    if movemean_length>0 %If smoothing and interpolation option was chosen by the user
        %Interpolate
        start=Track(k).time(1);
        stop=Track(k).time(end);
        t=start:dt:stop;
        tdoa=interp1(Track(k).time,Track(k).tdoa,t);

        %Smooth TDOAs
        tdoa_smooth=movmean(tdoa,movemean_length);
        tracks(k).time=t';
        tracks(k).time_local = t_serialdate(Track(k).ti(1):Track(k).ti(end));
        tracks(k).ti = Track(k).ti(1):Track(k).ti(end);
        tracks(k).tdoa=tdoa_smooth';
    else
        tracks(k).time=Track(k).time;
        tracks(k).time_local = t_serialdate(Track(k).ti);
        tracks(k).ti=Track(k).ti;
        tracks(k).tdoa=Track(k).tdoa;
    end


end

%remove tracks that are too short
indx=false(N,1);
for k=1:N
%     if length(tracks(k).time)>=min_tl
%         indx(k)=true;
%     end
    if (tracks(k).time(end)-tracks(k).time(1))>=min_tl
        indx(k)=true;
    end
end
tracks=tracks(indx);
end

