function [mb,wb,Pb,Tb,id] = birth_process_amplitude(model,Z,id)
%Create a GMM distribution based on measurements Z, then draw neborn states
%from that
%
%Inputs:
% - model = structure specifying birth intensity and covariance matrix
% - Z = measurements -a struct with two fields:
%   ~ measure.Z: a RFS of measurements- 1 x M cell array, where M is number  
%           of time steps. Each cell contains TDOA and amplitude of the 
%           cross-correlation information. 
%   ~ measure.T: a scalar specifying the number of time steps M 
% - id = a number to track max identity of the components (the highest 
%      label of the already existing components)
%
%Otuputs:
% - mb = states of the newborn targets - 2 x N matrix (N = number of
%   targets)
% - wb = weights of the newborn targets- 1 x N vector (N = number of
%   targets)
% - Pb = covariance matirces of the newborn targets - 2 x 2 x N array (N = 
%   number of targets)
% - Tb = identity of the newborn components - 1 x N vector (N = number of
%   targets)
% - id = a number specifying maximum identity of the components
%
%Pina Gruden, 2020

if ~isempty(Z)
    %set up mixture distribution for drawing Gaussian components:
    N=size(Z,2); %number of components in a Gaussian mixture
    Qb = repmat(0.1.*diag([model.Q(1),model.Q(4)]), [1,1,N]); %covariance matrix of the GMM components
    %With amlitude feature not all measurements are equally likely to be new targets
    %so take amplitudes to indicate likelihood of a measurement giving birth to
    %a new target (gmdistribution function: gm = gmdistribution(mu,sigma,p) -
    %and p is a mixing proportion - If p does not sum to 1, gmdistribution
    %normalizes it):
    Mub(:,1) = Z(1,:); %centeres (means) of the GMM components (based on TDOAS)
    switch model.bpv_type
        case 1
             Mub(:,2) = random(model.dottdoa_birth,N);
%                 Mub(:,2) = 0; %Uninformed prior
        case 2
            Zdeg=real(acosd((-model.c/model.d).*Z(1,:)));
            v=repmat(Zdeg',1,size(model.dottdoa_birth,2))-model.dottdoa_birth(1,:);
            [~,l]=min(abs(v),[],2);%index of the frequency that is closest to the measurement
            val=model.dottdoa_birth(2,l);
            Mub(:,2)=val;
    end

    objb = gmdistribution(Mub,Qb,Z(2,:));%mixing coefficients based on amplitudes- Z(2,:)
  
    Nb=100; %number of newborn components
    [mb,compInd] = random(objb,Nb); %states of the newborn targets
    mb=mb';
    sc = model.birth/sum(objb.ComponentProportion(compInd')); %get a scaler for multiplying the weights
    wb=sc.*(objb.ComponentProportion(compInd')); %weights of the newborn components (weighted proportional to amplitude)
    Pb = repmat([model.Q(1),0.5*model.Q(2);0.5*model.Q(3),model.Q(4)], [1,1,Nb]);  %covariance matirces of the newborn targets

    Tb=zeros(1,Nb);
    for j=1:Nb %number of Gaussian components in birth model
        Tb(j)=id;
        id=id+1;
    end
    
else
    mb=[];
    wb=[];
    Pb=[];
    Tb=[];
end
end

