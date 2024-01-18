%% Foot point from 2D origin O to uncertain line l
%
% x = sugr_construct_footpoint_Ol_Point_2D(l)
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% sw 9/2016
%
% See also sugr_Point_2D, sugr_get_isfinite_Point_2D, sugr_Line_2D,
% sugr_construct_mean_Point_2D, sugr_construct_midpoint_Point_2D,
% sugr_construct_intersection_Point_2D, sugr_construct_footpoint_xl_Point_2D

function Point_2D = sugr_construct_footpoint_Ol_Point_2D(l)

S3 = calc_S([0,0,1]');
A1 = calc_S(l.h) * S3;
h  = -A1 * l.h;
A = A1 + calc_S(S3 * l.h)';
Clhh = sugr_get_CovM_homogeneous_Vector(l);
Chh = A * Clhh * A';

% generate minimal parameters
Point_2D = sugr_Point_2D(h,Chh);
