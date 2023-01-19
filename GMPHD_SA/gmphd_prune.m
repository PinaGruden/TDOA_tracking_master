function [wn,mn,Pn,Tagn] = gmphd_prune(w,m,P,Tr,Tag)
% function gmphd_prune.m performs pruning of Gaussian components based on 
% a threshold Tr.
% 
% INPUTS
% - w - weights of the targets- 1 x N vector (N = number of targets)
% - m - states of the targets- d x N matrix (d = dimension of the state, N = number of targets)
% - P - covariance matirces of the targets - d x d x N array (d = dimension of the state,  
%   N = number of targets)
% - Tr - a number specifying the pruning threshold
% - Tag - identity of the components - 1 x N vector (N = number of
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


idx = find(w >= Tr); % Pruning- find targets with high enough weights

wn = w(idx).*(sum(w)/sum(w(idx))); %define new weights - Eq. (9.47), Mahler (2014), Advances in Statistical Multisource-mulititarget information fusion 
mn = m(:,idx);
Pn = P(:,:,idx);
Tagn=Tag(idx);
end

