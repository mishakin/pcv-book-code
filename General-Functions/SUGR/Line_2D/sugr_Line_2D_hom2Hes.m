%% Uncertain Hessian line from homogeneous line
%
% [e,Cee] = sugr_Line_2D_hom2Hes(l)
%
% * l    =  2D line struct {l.h,l.Crr}
%
% * e    = line parameters [phi, d]'
% * Cee  = 2 x 2 covariance matrix
%
%
% Jochen Meidow, Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% See also sugr_Line_2D, sugr_Line_2D_hom2cen,
% sugr_Line_2D_Hes2hom, sugr_Line_2D_Hes2cen, sugr_Line_2D_cen2hom, 
% sugr_get_Euclidean_Line_2D, sugr_get_centroid_Line_2D

function [e,Cee] = sugr_Line_2D_hom2Hes(ul)

l   = ul.h;
Cll = sugr_get_CovM_homogeneous_Vector(ul);

a = l(1);
b = l(2);
c = l(3);

p   = atan2(b,a);
s2  = a^2+b^2;
s   = sqrt(s2);
d   = -c/s;
e   = [p;d];

J   = [  -b/s2 ,   a/s2 ,   0  ;...
        c*a/s^3, c*b/s^3, -1/s];

Cee = J*Cll*J';

