function [lambda,Phi,omega] = DMD(X1,X2,eps,dt) 
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS: 
% X1 = X, data matrix
% X2 = X', shifted data matrix
% Columns of X1 and X2 are state snapshots 
% dt = time step advancing X1 to X2 (X to X')
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues

%% DMD
[U, S, V] = svd(X1, 'econ');
eval   = diag(S); 
[~,IX] = sort(eval,'descend');
IX     = IX(abs(eval)>eps); 

U_r = U(:, IX); % truncate to rank-r
S_r = S(IX, IX);
V_r = V(:, IX);
Atilde = U_r' * X2 * V_r *diag(1./diag(S_r)); % low-rank dynamics
[W_r, D] = eig(Atilde);
Phi = X2 * V_r *diag(1./diag(S_r)) * W_r*diag(1./diag(D)); % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

