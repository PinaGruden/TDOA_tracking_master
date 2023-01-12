function [m_predict,P_predict] = kalman_predict_multiple(model,m,P)      
% kalman_predict_multiple.m is a function that uses Kalman filter to
% predict target states according to system model.
%
% INPUTS:
% - model - a struct containing information on system model. 
%   Two fileds required: ~ F - system matrix  
%                        ~ Q - system noise covariance matrix
% - m - states of the targets- d x N matrix 
%   (d = dimension of the state, N = number of targets)
% - P - covariance matrices of the targets - d x d x N array 
%   (d = dimension of the state, N = number of targets)
%
% OUTPUTS:
% - m_predict - states of the predicted targets - d x N matrix 
%   (d = dimension of the state, N = number of targets)
% - P_predict - covariance matrices of the predicted targets - d x d x N 
%   array (d = dimension of the state, N = number of targets)
%   
%
%

plength= size(m,2);

m_predict = zeros(size(m));
P_predict = zeros(size(P));

for idxp=1:plength
    [m_temp,P_temp] = kalman_predict_single(model.F,model.Q,m(:,idxp),P(:,:,idxp));
    m_predict(:,idxp) = m_temp;
    P_predict(:,:,idxp) = P_temp;
end

function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)

m_predict = F*m;
P_predict =  Q+F*P*F'; 
