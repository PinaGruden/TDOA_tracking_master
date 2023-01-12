function [qz_temp,m_temp,P_temp] = kalman_update_single(z,H,R,m,P)
% kalman_update_single.m is a function that uses Kalman filter to
% update a single target state according to measurement model.
%
% INPUTS:
% - z - measurements for a given time step -  1 x N vector (N = number of
%   targets)
% - H - measurement matrix  
% - R - measurement noise covariance matrix
% - m - target state - d x 1 vector (d = dimension of the state)  
% - P - covariance matrix of the targets - d x d matrix
%   (d = dimension of the state)
%
% OUTPUTS:
% - qz_temp - target likelihood 1 x Z matrix (Z = number of measurements) 
% - m_temp - state of the updated target - d x Z matrix (d = dimension of 
%   the state,Z = number of measurements)
% - P_temp - covariance matrix of the updated target - d x d matrix
%   (d = dimension of the state)
%

mu = H*m;
S  = H*P*H'+R; 
Vs= chol(S); 
det_S= prod(diag(Vs))^2; 
inv_sqrt_S= inv(Vs); 
iS= inv_sqrt_S*inv_sqrt_S';
K  = P*H'*iS;

if size(z,1)==1
    qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*(z-repmat(mu,[1 size(z,2)])).*(iS*(z-repmat(mu,[1 size(z,2)]))));
else
    qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))));
end
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(mu,[1 size(z,2)]));
P_temp = (eye(size(P))-K*H)*P;
end