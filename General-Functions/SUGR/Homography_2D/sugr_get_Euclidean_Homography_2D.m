%% Euclideanly normalized 2D homography
% for quasi affine homographies w.r.t. disk around origin.
%
% [He,Cee] = sugr_get_isaffine_Homography(H)
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function [He,Cee] = sugr_get_Euclidean_Homography_2D(H)

if sugr_get_isquasiaffine_Homography_2D(H)
    He      = H.H/H.H(9);                            % Euclidean normalization
    Jhe     = 1/H.H(9) * [eye(8) -H.H(1:8)'/H.H(9)]; % Jacobian h->e
    Chh     = sugr_get_CovM_homogeneous_Homography_2D(H);
                                                     % CovM of vec H
    Cee     = Jhe * Chh * Jhe';                      % CovM of vec He
else
    He      = zeros(3);
    Cee     = zeros(9);
end
