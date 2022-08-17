function [tracks] = postprocess(Track,t_serialdate,parameters)

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
        tracks(k).tdoa=tdoa_smooth';
    else
        tracks(k).time=Track(k).time;
        tracks(k).time_local = t_serialdate(Track(k).ti);
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

