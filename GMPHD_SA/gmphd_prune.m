function [wn,mn,Pn,Tagn] = gmphd_prune(w,m,P,Tr,Tag)
%Pruning of Gaussian components based on a threshold Tr

idx = find(w >= Tr); % Pruning- find targets with high enough weights

wn = w(idx).*(sum(w)/sum(w(idx))); %define new weights - Eq. (9.47), Mahler (2014), Advances in Statistical Multisource-mulititarget information fusion 
mn = m(:,idx);
Pn = P(:,:,idx);
Tagn=Tag(idx);
end

