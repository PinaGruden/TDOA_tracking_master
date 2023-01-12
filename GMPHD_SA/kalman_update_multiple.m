function [qz_update,m_update,P_update,Tag_update] = kalman_update_multiple(z,model,m,P,Tag)
% kalman_update_multiple.m is a function that uses Kalman filter to
% update target states according to measurement model.
%
% INPUTS:
% - z - measurements for a given time step -  1 x N vector (N = number of
%   targets)
% - model - a struct containing information on system model. 
%   Two fileds required: ~ H - measurement matrix  
%                        ~ R - measurement noise covariance matrix
% - m - states of the targets- d x N matrix 
%   (d = dimension of the state, N = number of targets)
% - P - covariance matrices of the targets - d x d x N array 
%   (d = dimension of the state, N = number of targets)
% - Tag - identities (labels) of the targets - 1 x N vector (N = number of
%   targets)
%
% OUTPUTS:
% - qz_update - multitarget likelihood N x Z matrix (N=number of targets, Z
%   = number of measurements)
% - m_update - states of the updated targets - d x N x Z matrix 
%   (d = dimension of the state, N = number of targets, Z = number of 
%   measurements)
% - P_update - covariance matrices of the updated targets - d x d x N 
%   array (d = dimension of the state, N = number of targets)
% - Tag_update - identities (labels) of the updated targets - Z x N vector 
%   (Z = number of measurements, N = number of targets)

plength= size(m,2);
zlength= size(z,2);
x_dim=size(m,1);


qz_update= zeros(plength,zlength);
m_update = zeros(x_dim,plength,zlength);
P_update = zeros(x_dim,x_dim,plength);

for id=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single(z,model.H,model.R,m(:,id),P(:,:,id));
        qz_update(id,:)   = qz_temp;
        m_update(:,id,:) = m_temp;
        P_update(:,:,id) = P_temp;        
end

Tag_update = repmat(Tag,zlength,1);

end



