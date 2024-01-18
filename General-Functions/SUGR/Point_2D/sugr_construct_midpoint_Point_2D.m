%% Mid point of two 2D points
%
% m = sugr_construct_midpoint_Point_2D(x,y)
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
% sw 9/2016
%
% See also sugr_Point_2D, sugr_get_isfinite_Point_2D, sugr_Line_2D,
% sugr_construct_mean_Point_2D, sugr_construct_intersection_Point_2D,
% sugr_construct_footpoint_Ol_Point_2D, sugr_construct_footpoint_xl_Point_2D

function Point_2D = sugr_construct_midpoint_Point_2D(x,y)

A = calc_CentralMatrix(y.h);
B = calc_CentralMatrix(x.h);
h   = A * x.h;
% hc  = B * y.h;
Cxhh = sugr_get_CovM_homogeneous_Vector(x);
Cyhh = sugr_get_CovM_homogeneous_Vector(y);
Chh =  A * Cxhh * A' + B * Cyhh * B';

% generate reduced parameters
Point_2D = sugr_Point_2D(h,Chh);
