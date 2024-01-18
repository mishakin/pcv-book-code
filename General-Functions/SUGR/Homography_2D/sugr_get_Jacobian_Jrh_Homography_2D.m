%% Get Jacobian Jrh dr/dh for estimating 2D Homography  
%
% [Jrh,Jsh,Jks,Jrk] = sugr_get_Jacobian_Jhr_Homography_2D(H)
%
% H      3x3-matrix spectrally normalized homography
%
% Jrh    8x9-matrix Jacobian dr/dh, h -> r
% Jsh    9x9-matrix Jacobian ds/sh, h -> s
% Jks    9x9-matrix Jacobian dk/ds, s -> k
% Jrk    8x9-matrix Jacobian dr/dk, k -> r
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 3/2011, 
% see PCV pp. 385 ff. (r == Delta p)
% and App. A.15

function [Jrh,Jsh,Jks,Jrk] = sugr_get_Jacobian_Jrh_Homography_2D(H)


Hit    = inv(H');                      % H^-T for Jacobian   
vHit   = Hit(:);                       % vectorize
Jsh    = eye(9) - 1/3 * H(:)*vHit';    % Jacobian J: ds/dh, h -> s
Jks    = kron(Hit,eye(3));             % Jacobian J: dk/ds, s->k (Hit otimes I_3)    
Jrk    = [eye(8),zeros(8,1)];          % Jacobian J: dr/dk, k->h

Jrh    = Jrk * Jks * Jsh;              % Jacobian J: dr/dr, h->r
