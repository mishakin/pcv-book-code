% sugr_construct_mean_Point_2D: mid point of two 2D points
%
% m = sugr_construct_mean_Point_2D(x,y)
%
% Wolfgang Förstner 7/2014
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_get_isfinite_Point_2D, sugr_Line_2D,
% sugr_construct_midpoint_Point_2D, sugr_construct_intersection_Point_2D,
% sugr_construct_footpoint_Ol_Point_2D, sugr_construct_footpoint_xl_Point_2D

function Point_2D = sugr_construct_mean_Point_2D(x,y)

% use mid point as approximate value
xa = sugr_construct_midpoint_Point_2D(x,y);

% use tangent space
J = null(xa.h');

% reduced coordinates
xr = J'*x.h;
yr = J'*y.h;
Jx = J'*calc_Rot_ab(x.h,xa.h)*null(x.h');
Jy = J'*calc_Rot_ab(y.h,xa.h)*null(y.h');
Wxx = inv(Jx*x.Crr*Jx');
Wyy = inv(Jy*y.Crr*Jy');

% mean in tangent space
Cmrmr = inv(Wxx+Wyy);
xmr   = Cmrmr*(Wxx*xr+Wyy*yr);                                             %#ok<*MINV>

% go to homogeneous coordinates
Cmm = J * Cmrmr * J';                            
xm = xa.h + J*xmr;
xm = xm/norm(xm); 

% generate reduced parameters
Point_2D = sugr_Point_2D(xm,Cmm);
