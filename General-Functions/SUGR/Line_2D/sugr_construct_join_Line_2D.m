%% Construct 2D line joining two points x and y 
%
% Line_2D = sugr_construct_join_Line_2D(x,y)
%
% * x,y =  2D points, struct {x.h,x.Crr}
% * Line_2D =  2D line, struct {l.h,l.Crr}

% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D, sugr_construct_parallel_lx_Line_2D, 
% sugr_construct_parallel_lO_Line_2D, sugr_construct_midline_Line_2D, 
% sugr_construct_footpoint_Ol_Point_2D, sugr_construct_footpoint_xl_Point_2D, 
% sugr_construct_join_Line_3D


function Line_2D = sugr_construct_join_Line_2D(x,y)
     
        Sx =   calc_S(x.h);
        Sy = - calc_S(y.h);
        h   = Sx * y.h;
        Cxhh = sugr_get_CovM_homogeneous_Vector(x);
        Cyhh = sugr_get_CovM_homogeneous_Vector(y);
        Chh = Sx * Cyhh * Sx' + Sy * Cxhh * Sy'; 
        
% generate line
Line_2D = sugr_Line_2D(h,Chh);
