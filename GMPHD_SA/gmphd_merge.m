function [wn,mn,Pn,Tagn] = gmphd_merge(w,m,P,U,Tag)
% function gmphd_merge.m merges Gaussian components that are close together
% based on the threshold U.
%
% INPUTS:
% - w - weights of the targets- 1 x N vector (N = number of targets)
% - m - states of the targets- d x N matrix (d = dimension of the state, N = number of targets)
% - P - covariance matirces of the targets - d x d x N array (d = dimension of the state,  
%   N = number of targets)
% - U - a number specifying the merging threshold
% - Tag - identity of the components - 1 x N vector (N = number of
%   targets)
%
% OUTPUTS:
% - wn - weights of the merged targets
% - mn - states of the merged targets
% - Pn - covariance matirces of the merged targets
% - Tagn - identities of the merged components
%
%
% Pina Gruden


xdim=size(m,1);
In=1:length(w);
l=0; 
%preallocate (to the max possible size, delete empty spaces later)
wn=zeros(1,length(In)); mn=zeros(xdim,length(In));Pn=zeros(xdim,xdim,length(In));
% Tagn=cell(1,length(In));
Tagn=zeros(1,length(In));

while ~isempty(In)
    l=l+1;
    [~,j] = max(w(In)); %find biggest of the pruned weights
    %now find indices of the means of Gaussian components that are within
    %distance U from the Gaussian component with highest weight
    %Use Mahalanobis distance
    difs=m(:,In)-repmat(m(:,In(j)),size(In));
    mahal=diag(difs'*(pinv(P(:,:,In(j)))*difs))'; %replaced (P(:,:,In(j))\difs with pinv to avoid ill conditioning
%     Pinv=inv(P(:,:,In(j)));
%     mahal=zeros(1,length(In));
%     for n=1:length(In)
%         mahal(n)= (m(:,In(n))-m(:,In(j)))'*Pinv*(m(:,In(n))-m(:,In(j)));
%     end
    L=  In(mahal<=U);%get indices
    
    %calculate new weights,means and covariances for merged components
    wn(l)=sum(w(L));
    mn(:,l)=sum(w(L).*m(:,L),2)/wn(l); 
    P_l=zeros(xdim,xdim);
    for i=1:length(L)
        P_l = P_l + w(L(i))*(P(:,:,L(i))+(mn(:,l) - m(:,L(i)))*(mn(:,l) - m(:,L(i)))');
    end
    Pn(:,:,l) = P_l;
    Pn(:,:,l)= Pn(:,:,l)/wn(l);
    
%     Tagn{l}=Tag{In(j)};
    Tagn(l)=Tag(In(j));
    
    % remove elements of L from I (I:= I\L)
    matches=zeros(size(L));
    for ix=1:length(L),matches(ix)=find(L(ix)==In);end
    In(matches)=[];
end
%now delete empty spaces
wn=wn(1:l);mn=mn(:,1:l);Pn=Pn(:,:,1:l);
Tagn=Tagn(1:l);

end

