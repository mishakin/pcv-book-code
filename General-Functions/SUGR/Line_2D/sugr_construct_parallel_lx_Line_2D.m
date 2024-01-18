%% Construct 2D line parallel to line l through point x
%
% Line_2D = sugr_construct_parallel_lx_Line_2D(l,x)
%
% * l =  2D line struct {l.h,l.Crr}
% * x =  2D point struct {x.h,x.Crr}
%
% * Line_2D  =  2D line struct {l.h,l.Crr}
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_construct_parallel_lO_Line_2D,
% sugr_construct_midline_Line_2D, sugr_construct_footpoint_Ol_Point_2D,
% sugr_construct_footpoint_xl_Point_2D, 
% sugr_construct_join_Line_2D, sugr_construct_join_Line_3D

function Line_2D = sugr_construct_parallel_lx_Line_2D(l,x)

        S3 = calc_S([0,0,1]');
        A = -calc_S(x.h)*S3;
        B =  calc_S(S3*l.h);
        h   = A * l.h;
        Cxhh = sugr_get_CovM_homogeneous_Vector(x);
        Clhh = sugr_get_CovM_homogeneous_Vector(l);
        Chh = A * Clhh * A' + B * Cxhh * B'; 
     
% generate line
Line_2D = sugr_Line_2D(h,Chh);
