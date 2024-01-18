%% Construct line parallel to line l through origin
%
% lp = sugr_construct_parallel_lO_Line_2D(l)
%
% * l, lp =  2D line struct {l.h,l.Crr}
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_construct_parallel_lx_Line_2D, 
% sugr_construct_midline_Line_2D, sugr_construct_footpoint_Ol_Point_2D,
% sugr_construct_footpoint_xl_Point_2D, 
% sugr_construct_join_Line_2D, sugr_construct_join_Line_3D


function Line_2D = sugr_construct_parallel_lO_Line_2D(l)
   
        G3 = diag([1,1,0]);
        h   = G3 * l.h;
        Clhh = sugr_get_CovM_homogeneous_Vector(l);
        Chh = G3 * Clhh * G3; 
    
% generate line
Line_2D = sugr_Line_2D(h,Chh);

