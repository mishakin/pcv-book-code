%% Test whether homography is affine
% preimage of line at infinity: H(3,1:2)=[0,0]
%
% sugr_get_isaffine_Homography_2D: boolean = transformation is affinity
%
% isaffine = sugr_get_isaffine_Homography_2D(H)
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function isaffine = sugr_get_isaffine_Homography_2D(H)

global Threshold_Euclidean_Normalization

isaffine = norm(H.H(3,1:2)) < Threshold_Euclidean_Normalization;
