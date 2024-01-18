%% Show information of 3D point on console
%
% sugr_show_Point_3D(Point_3D, [name])

% * Point_3D = sugr object, structure 
%      .h     = spherically normalized homogeneous coordinates
%      .Crr   = reduced covariance matrix
%      .Jr    = null space of .h'
% * name = Name of object, optional
%
% Wolfgang Förstner 1/2011
% wfoerstn@uni-bonn.de
%
% sw 9/2016
% sw 1/2018
%
% See also sugr_Point_3D

function sugr_show_Point_3D(x, name)

if nargin<2
    name = 'X';
end

global type_name

f    = sugr_get_isfinite_Point_3D(x);
typ_no   = x.type;

if f
    fprintf('\n%s: finite %s\n',name,type_name{typ_no})
else
    fprintf('\n%s: infinite %s\n',name,type_name{typ_no})
end

[e,Cee] = sugr_get_Euclidean_Point_3D(x);
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\n',             e(1),Cee(1,1), Cee(1,2), Cee(1,3))
fprintf('\t%s_e =\t%5.3f\t\tCov_ee =\t%5.3f %5.3f %5.3f\n',name,e(2),Cee(2,1), Cee(2,2), Cee(2,3))
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\n',             e(3),Cee(3,1), Cee(3,2), Cee(3,3))

h = x.h;
Chh  = sugr_get_CovM_homogeneous_Vector(x);
Crr  = x.Crr;
fprintf('\n\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f %5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\n',                     h(1),Chh(1,1), Chh(1,2), Chh(1,3), Chh(1,4), Crr(1,1), Crr(1,2), Crr(1,3))
fprintf('\t%s_h =\t%5.3f\t\tCov_hh =\t%5.3f %5.3f %5.3f %5.3f\t\tCov_rr =\t%5.3f %5.3f %5.3f\n',name,h(2),Chh(2,1), Chh(2,2), Chh(2,3), Chh(2,4), Crr(2,1), Crr(2,2), Crr(2,3))
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f %5.3f\t\t\t\t\t%5.3f %5.3f %5.3f\n',            h(3),Chh(3,1), Chh(3,2), Chh(3,3), Chh(3,4), Crr(3,1), Crr(3,2), Crr(3,3))
fprintf('\t\t\t%5.3f\t\t\t\t\t%5.3f %5.3f %5.3f %5.3f\n',                                      h(4),Chh(4,1), Chh(4,2), Chh(4,3), Chh(4,4))



