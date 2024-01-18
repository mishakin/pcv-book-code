%% Test whether homography is quasi affine
% boolean = oriented in disk of radius
% preimage of line at infinity must have distance larger radius
%
% isquasiaffine = sugr_get_isquasiaffine_Homography(H)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function isquasiaffine = sugr_get_isquasiaffine_Homography_2D(H)


global Threshold_preimage_line_infinity

radius = Threshold_preimage_line_infinity;

isquasiaffine = norm(H.H(3,1:2)) < radius * abs(H.H(3,3));

