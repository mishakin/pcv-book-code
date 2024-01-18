%% Construct midline of two 2D points x and y
%
% Line_2D = sugr_construct_midline_Line_2D(x,y)

% * x,y =  2D points, struct {x.h,x.Crr}
% * Line_2D =  2D line, struct {l.h,l.Crr}
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_2D, sugr_construct_parallel_lx_Line_2D,  
% sugr_construct_parallel_lO_Line_2D, sugr_construct_footpoint_Ol_Point_2D, 
% sugr_construct_footpoint_xl_Point_2D, sugr_construct_join_Line_2D, 
% sugr_construct_join_Line_3D

function Line_2D = sugr_construct_midline_Line_2D(x,y)

G3 = diag([1,1,0]);
z = calc_CentralMatrix(x.h)*y.h;
l = calc_S(x.h)*y.h;
A = -calc_S(G3 * l) * calc_CentralMatrix(y.h) ...
    -calc_S(z) * G3 * calc_S(y.h);
B = -calc_S(G3 * l) * calc_CentralMatrix(x.h) ...
    +calc_S(z) * G3 * calc_S(x.h);
h = +calc_S(z) * G3 * calc_S(x.h) * y.h;
Cxhh = sugr_get_CovM_homogeneous_Vector(x);
Cyhh = sugr_get_CovM_homogeneous_Vector(y);
Chh = A * Cxhh * A' + B * Cyhh * B';


% generate line
Line_2D = sugr_Line_2D(h,Chh);
