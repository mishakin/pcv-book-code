%% Centroid representation from minimal representation
%
% [x0,p,sp,sq] = sugr_get_centroid_Line_2D(l)
%
% * l     = 2D line struct, minimal representation {l.h,l.Crr}
%
% * x0 = centroid, [x-coord. y-coord.]'
% * p  = direction of normal
% * sp = standard deviation of normal direction
% * sq = standard deviation across line
%
% Jochen Meidow, Wolfgang Förstner
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [x0,p,sp,sq] = sugr_get_centroid_Line_2D(l)

[e,Cee]      = sugr_get_Euclidean_Line_2D(l);
[x0,p,sp,sq] = sugr_Line_2D_Hes2cen(e,Cee);
