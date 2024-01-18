%% sugr_show_Homography_2D(H)
%
% show elements of homography
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 1/2011

function sugr_show_Homography_2D(H)

[He,~] = sugr_get_Euclidean_Homography_2D(H);

He    = He                                                                
detH  = det(H.H)
isquasiaffine = sugr_get_isquasiaffine_Homography_2D(H)
isaffine      = sugr_get_isaffine_Homography_2D(H)
chi           = sugr_get_chirality_Homography_2D(H)
