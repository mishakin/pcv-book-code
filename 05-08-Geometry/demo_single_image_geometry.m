%% DEMO single image geometry
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 07/17
% wenzel@igg.uni-bonn.de

clc
close all

addpath(genpath('../General-Functions'))

fprintf('\n------ DEMO Single Image Geometry ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')

%% Given

fprintf('\nGiven:\n')

% object point
fprintf('\n3D object point\n')
X = [3,2,5,6]'; disp(['X = [', num2str(X'./X(end)),']'])

% Interior orientation, calibration parameters
c = 1;  % camera constant
s = 0;  % shear
xH = 0; % principal point
yH = 0;
m = 0;  % scale difference

% exterior orentation
fprintf('\nProjection center\n')
X0 = [1,2,3,1]'; disp(['XO = [', num2str(X0'),']'])
% rotation, rotation vector (r tan(phi/2))
r=[0,0,0]';


fprintf('\nInterior orientaion, calibration matrix\n')
K = [c,s*c,xH;0,c*(1+m),yH;0,0,1]                                          %#ok<NOPTS>

fprintf('\nExterior orientaion, rotation matrix\n')
R = calc_Rot_r(r)                                                          %#ok<NOPTS>

%% Projection matrix for points, P = K R [I | -XO]
fprintf('\nProjection matrix for points, P = K R [I | -XO]\n')
P = calc_P_from_KRZ(K,R,X0(1:3)/X0(4))                                     %#ok<NOPTS>

%% Projection matrix for 3D-lines, see PCV Eq.(12.71)
fprintf('\nProjection matrix for 3D-lines, see PCV Eq.(12.71)\n')
Q = calc_Q_from_P(P)                                                       %#ok<NOPTS>

%% Projection of image point
fprintf('\nProjection of object point --> homogeneous coordiantes of image point\n')
xh = P*X; disp(['x_h = [', num2str(xh'),']'])

fprintf('\nEuclidian coordinates of image point\n')
xe = xh(1:2)/xh(3); disp(['x = [', num2str(xe'),']'])

%% Projection ray in Plücker coordiantes
fprintf('\nProjection ray in Plücker coordiantes\n')
Lxs = calc_Dual*Q'*xh ; disp(['L_xs = [', num2str(Lxs'),']'])

fprintf('\nCheck  whether projection ray passes through 3D-point: should be 0\n')
test_val = calc_Gammadual(Lxs) * X ; disp(['[', num2str(test_val'),']'])


