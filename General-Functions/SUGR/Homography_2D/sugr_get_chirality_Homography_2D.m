%% etermines chirality of quasi affine homography (within disk)
% chirality of quasi affinity +-1 or 0 = undefined
%
% chi = sugr_get_chirality_Homography_2D(H,radius)
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2011

function chi = sugr_get_chirality_Homography_2D(H)

if sugr_get_isquasiaffine_Homography_2D(H)
            chi = sign(H.H(3,3)); 
else
    chi = 0;
end
