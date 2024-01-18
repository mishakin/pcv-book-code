% demo_Point_3D: test routine for generating 3D point
%
% Point_3D = structure 
%
% *           .h = spherically normalized homogeneous coordinates
% *           .Crr = half vector of reduced covariance matrix
% *           .type = 'Point_2D'
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_3D

clearvars
clc

%% initialize sugr
% 
addpath(genpath('../General-Functions/'));
sugr_INIT

%% Sorted according to number of arguments

%% (1) Euclidean point vector without CovM (==0)
Xe = [1,2,3]';
X = sugr_Point_3D(Xe);
sugr_show_Point_3D(X,'X');

%% (1) Homogeneous point vector without CovM (==0)
Xh = [1,2,3,2]';
X = sugr_Point_3D(Xh);
sugr_show_Point_3D(X,'X');

%% (2) Homogeneous point vector with CovM
Xh = [1,2,3,2]';
Chh = [3,2,1,0;2,3,2,1;1,2,3,2;0,1,2,3]*0.1;
X = sugr_Point_3D(Xh,Chh);
sugr_show_Point_3D(X,'X');

%% (2) Homogeneous vanishing point vector
Xh = [3,4,2,0]';
Chh = [3,2,1,0;2,3,2,1;1,2,3,2;0,1,2,3]*0.1;
X = sugr_Point_3D(Xh,Chh);
sugr_show_Point_3D(X,'X');

%% (6) Euclidean point coordinates with with standard dev
x = 1;
y = 2;
z = 3;
sx = 0.1;
sy = 0.3;
sz = 0.2;
X=sugr_Point_3D(x,y,z,sx,sy,sz);
sugr_show_Point_3D(X,'X');

