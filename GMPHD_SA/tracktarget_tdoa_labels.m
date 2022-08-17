function [ Track ] = tracktarget_tdoa_labels(XTag,Xk,model)
%tracking labels is a function that finds targets with the same label and
%collects them into tracks.

%Pina Gruden, 2020

%based on the tracking_labels_v2.m from my PhD (from the SMCPHD)

IDlist= unique([XTag{:}]);
N=numel(IDlist); %number of all tracks
if N==0 %there are no tracks
    Track=[];
    return
end
%preallocate
Track(N).tdoa =[];
Track(N).dottdoa =[];
Track(N).time =[];
Track(N).label =[];
Track(N).ti =[];

t=0:model.dt:size(XTag,1)*model.dt-model.dt;

for k=1:size(XTag,1)

    if ~isempty(XTag{k})
        [ni, indxs] = histc(XTag{k}, IDlist); %gives indices of which IDs the tags belong to
        %Do NOT use histcounts- [ni2,~,indxs2] = histcounts(XTag{k},
        %IDlist)- because the function bins together the last two nins
        %hence if the label corresponds to the last number on the IDlist,
        %it will be given a wrong index (indxs) that will correspond to the
        %IDlist(end)-1.
        multiple = find(ni > 1); %gives indication wether there are estimates with the same ID
        
        %identify which IDs are repeated and which not
        ID=ismember(indxs, multiple); %Logical array
        
        %~~~~~~~~~~~~~~ For IDS that are NOT repeated  ~~~~~~~~~~~~~~~~~~
        indx=indxs(ID==0);
        for m=1:length(indx)
            Track(indx(m)).tdoa=[Track(indx(m)).tdoa;Xk{k}(1,indxs==indx(m))];
            Track(indx(m)).dottdoa=[Track(indx(m)).dottdoa;Xk{k}(2,indxs==indx(m))];
            Track(indx(m)).time=[Track(indx(m)).time; t(k)];
            Track(indx(m)).label=[Track(indx(m)).label;IDlist(indx(m))];
            Track(indx(m)).ti = [Track(indx(m)).ti,k]; %time step index
        end
        
        %~~~~~~~~~~~~~~ For IDS that ARE repeated  ~~~~~~~~~~~~~~~~~~
        if ~isempty(multiple)
            for m=1:length(multiple)
            %------------ Resolve based on distance ------------------
            
            detections = Xk{k}(:,indxs==multiple(m));
            %predictions:
            M=multiple(m);
            
            if isempty(Track(M).ti) 
              %if this is the first time this track/ID is detected, then no 
              %info is there,so just take a mean between the estimates to 
              %be the real estimate  
              meanest=mean(detections,2);
              Track(M).tdoa=[Track(M).tdoa;meanest(1)];
              Track(M).dottdoa=[Track(M).dottdoa;meanest(2)];
              Track(M).time=[Track(M).time; t(k)];
              Track(M).label=[Track(M).label;IDlist(M)];
              Track(M).ti = [Track(M).ti,k];
            else
                
                %compute tdoa derivative : using tdoa derivative
                if  numel(Track(M).ti) <= 5
                    alpha =median(Track(M).dottdoa);
                else
                    alpha = median(Track(M).dottdoa(end-5:end));
                end
                
                %deterime how many time steps ago it was last detected
                Nstp=k-Track(M).ti(end);
                A=model.F; A(1,2)=A(1,2)*Nstp;
                prediction=A*[Track(M).tdoa(end);alpha];
                
                dif=detections-prediction;
                %compute Mahalanobis distance
                sig = [model.Q(1),0;0,model.Q(4)];
                mah = sqrt(sum(dif'.^2/sig,2));
                [~,i]=min(mah); %index of the component with min distance from expected location
                
                Track(M).tdoa=[Track(M).tdoa;detections(1,i)];
                Track(M).dottdoa=[Track(M).dottdoa;detections(2,i)];
                Track(M).time=[Track(M).time; t(k)];
                Track(M).label=[Track(M).label;IDlist(M)];
                Track(M).ti = [Track(M).ti,k]; %time step index
            end
            end
        end
    end
end
end



