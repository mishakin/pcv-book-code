%% Foot point from Point_2D x onto line l
%
% x = sugr_construct_footpoint_xl_Point_2D(x,l)
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
% 
% sw 9/2016
%
% See also sugr_Point_2D, sugr_get_isfinite_Point_2D, sugr_Line_2D,
% sugr_construct_mean_Point_2D, sugr_construct_midpoint_Point_2D,
% sugr_construct_intersection_Point_2D,sugr_construct_footpoint_Ol_Point_2D,

function Point_2D = sugr_construct_footpoint_xl_Point_2D(x,l)

G3 = diag([1,1,0]);
Sl = calc_S(l.h);
Sx = calc_S(x.h);
A1 = Sl * Sx * G3;
h   = A1 * l.h;
A =   A1 + calc_S(Sx * G3 * l.h)';
B =  -Sl * calc_S(G3 * l.h);
Clhh = sugr_get_CovM_homogeneous_Vector(l);
Cxhh = sugr_get_CovM_homogeneous_Vector(x);
Chh = A * Clhh * A' + B * Cxhh * B';

%generate reduced parameters
Point_2D = sugr_minimal_vector(h,Chh);
