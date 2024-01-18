%% Uncertain centered line from uncertain homogeneous line
% 
% [x0,p,sp,sq] = sugr_Line_2D_hom2cen(l)
%
% * l    =  2D line struct {l.h,l.Crr}
%
% * x0 = centroid, [x-coord. y-coord.]'
% * p  = direction of normal
% * sp = standard deviation of normal direction
% * sq = standard deviation across line
%
% Jochen Meidow, Wolfgang Förstner 
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_Line_2D_hom2Hes, 
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [x0,p,sp,sq] = sugr_Line_2D_hom2cen(l)

[e,Cee]      = sugr_Line_2D_hom2Hes(l);
[x0,p,sp,sq] = sugr_Line_2D_Hes2cen(e,Cee);

