%% Show/print all representations of 2D line on console
%
% sugr_show_Line_2D(Line_2D, [name])
%
% * Line_2D = sugr object, structure 
%      .h     = spherically normalized homogeneous coordinates
%      .Crr   = reduced covariance matrix
%      .Jr    = null space of .h'
% * name = Name of object, optional
%
% Wolfgang Förstner  1/2011
% wfoerstn@uni-bonn.de
%
% sw 1/2018
%
% See also sugr_show_Point_2D, sugr_show_Point_3D, sugr_show_Homography_2D

function sugr_show_Line_2D(Line_2D,name)

if nargin<2
    name = 'l';
end

global type_name

f    = sugr_get_isfinite_Line_2D(Line_2D);
typ_no   = Line_2D.type;

if f
    fprintf('\n%s: finite %s\n',name,type_name{typ_no})
else
    fprintf('\n%s: infinite %s\n',name,type_name{typ_no})
end


[e,Cee] = sugr_get_Euclidean_Line_2D(Line_2D);
fprintf('\t%s_e =\t%5.3f\t\tCov_ee =\t%5.3f %5.3f\n',name,e(1),Cee(1,1), Cee(1,2))
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f\n',e(2),Cee(2,1), Cee(2,2))

h = Line_2D.h;
Chh  = sugr_get_CovM_homogeneous_Vector(Line_2D);
Crr  = Line_2D.Crr;
fprintf('\n\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\n',h(1),Chh(1,1), Chh(1,2), Chh(1,3))
fprintf('\t%s_h =\t%5.3f\t\tCov_hh =\t%5.3f %5.3f %5.3f\t\tCov_rr =\t%5.3f %5.3f\n',name,h(2),Chh(2,1), Chh(2,2), Chh(2,3), Crr(1,1), Crr(1,2))
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\t\t\t\t\t%5.3f %5.3f\n',h(3),Chh(3,1), Chh(3,2), Chh(3,3), Crr(2,1), Crr(2,2))

[x0,phi,sp,sq] = sugr_get_centroid_Line_2D(Line_2D);
fprintf('\n\tx_0 = [%5.3f, %5.3f], phi = %5.3f°, sp = %5.3f°, sq = %5.3f\n',x0,phi*180/pi,sp*180/pi,sq)



