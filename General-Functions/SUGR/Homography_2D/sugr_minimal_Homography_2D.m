%% minimal representation of homography
% i.e., spectrally normalize uncertain homography
%
% [H,Jhr] = sugr_minimal_Homography_2D(Ha,Chh)
%
% * Ha     uncertain regular homography
% * Chh    covariance matrix of vec Ha
%
% H      minimal representation:
% * .H     homography matrix, sprectrally normalized
% * .Crr   covariance matrix of vector of left differential factor
%
% Jrh      Jacobian dr/dh, h -> r
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 3/2011

function [H,Jrh] = sugr_minimal_Homography_2D(Ha,Chh)

det_H  = abs(det(Ha));              % positiv determinant
factor = det_H^(-1/3);              % normalizing factor
H.H    = Ha  * factor;              % normalization: det = +-1
Chh    = Chh * factor^2;

% Jhr    % Jacobian J: h->r
Jrh = sugr_get_Jacobian_Jrh_Homography_2D(H.H);

H.Crr  = Jrh * Chh * Jrh';           % CovM of reduced vector
H.type = 20;       
