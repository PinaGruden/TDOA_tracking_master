function [wn,mn,Pn,Tagn] = gmphd_limit(w,m,P,Jmax,Tag)
% function gmphd_limit.m eliminates Gaussian components based on their
% weights w, so that the final number does not exceed Jmax components
%
% INPUTS
% - w - weights of the targets- 1 x N vector (N = number of targets)
% - m - states of the targets- d x N matrix (d = dimension of the state, N = number of targets)
% - P - covariance matirces of the targets - d x d x N array (d = dimension of the state,  
%   N = number of targets)
% - Jmax - number specifying the maximum allowed number of Gaussian
%   components
% - Tag - identities of the components - 1 x N vector (N = number of
%   targets)
%
% OUTPUTS:
% - wn - weights of the pruned targets
% - mn - states of the pruned targets
% - Pn - covariance matirces of the pruned targets
% - Tagn - identities of the pruned components
%
%
% Pina Gruden


if numel(w) > Jmax
    [~,ranking] = sort(w,'descend'); %gives indices of ranking in descending order
    indx=ranking(1:Jmax); %indices that we want to keep 
    wn = w(indx); 
    wn=wn*(sum(w)/sum(wn));%re-weight the components
    mn=m(:,indx);
    Pn= P(:,:,indx);
    Tagn=Tag(indx);
else
    wn=w;
    mn=m;
    Pn=P;
    Tagn=Tag;
end

end

