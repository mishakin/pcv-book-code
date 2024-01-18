%% determines constraint values and Jacobians for 2D similarity
%
% Model:
% cg(x,l) = zr - s*zl - t = 0
%
% [l(3)]    - [x(1) -x(2)] [l(1)] - [x(3)] = [0]
% [l(4)]    - [x(2)  x(1)] [l(2)] - [x(4)]   [0]
%
% l         = 4x1-vector [xl yl,xr yr] = [zl, zr]
% x         = 4x1-vector [ a  b, c  d] = [ s,  t]
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de

function [cg,At,Bt] = cg_similarity_2D(l,la,xa)

At  = [-la(1)  la(2) -1  0 ;...
    -la(2) -la(1)  0 -1];
Bt  = [-xa(1)  xa(2) 1 0;...
    -xa(2) -xa(1) 0 1];

cg = -( [la(3);...
    la(4)] - ...
    [xa(1) -xa(2);
    xa(2)  xa(1)] * ...
    [la(1);...
    la(2)] - ...
    [xa(3); ...
    xa(4)])...
    + Bt * (la-l);
