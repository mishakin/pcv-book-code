%% Intersection 2D point two lines l and m
%
% x = sugr_construct_interection_Point_2D(l,m)
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
% 
% sw 9/2016
%
% See also sugr_Point_2D, sugr_get_isfinite_Point_2D, sugr_Line_2D,
% sugr_construct_mean_Point_2D, sugr_construct_midpoint_Point_2D,
% sugr_construct_footpoint_Ol_Point_2D, sugr_construct_footpoint_xl_Point_2D

function Point_2D = sugr_construct_intersection_Point_2D(l,m)

Sl = calc_S(l.h);
Sm = calc_S(m.h);
h   = Sl * m.h;
Clhh = sugr_get_CovM_homogeneous_Vector(l);
Cmhh = sugr_get_CovM_homogeneous_Vector(m);
Chh = Sl * Cmhh * Sl' + Sm * Clhh * Sm';

% generate reduced parameters
Point_2D = sugr_Point_2D(h,Chh);
