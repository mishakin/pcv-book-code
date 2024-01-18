%% test mirror image
% example according to Figure 13.11 page 579
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

addpath(genpath('../General-Functions'))

close all
clc

% three points: X=right, Y=depth, Z=heigth
disp('----------- mirror image -----------')
N = 3;

%% in original position, will be rotated later
% points
disp('Create points and plane in generic situation ...')
disp('original points')
X0 = [1.0,1.0,1.2;
      1.2,1.0,1.5;
      1.1,3.0,2.0]                                                         %#ok<*NOPTS>

% mirror plane
disp('mirrow plane')
A0 = [1,0,0,-2]'

% mirrored points
disp('mirrowed points')
Y0 = [4-X0(:,1),X0(:,2),X0(:,3)]

%% now rotate to obtain generic situation
% rotation
disp('rotation to obtain generic situation')
Rp = calc_Rot_rod([0,0,0.05]')                                             

% rotated points, plane
disp('rotated points (see Tab. 13.4)')
X = X0*Rp'
Y = Y0*Rp'
disp('transformation matrix for plane')
M = [Rp, zeros(3,1);zeros(1,3),1]
Ma = adjunctMatrix(M);
disp('rotated mirrow plane')
A = Ma' * A0

% check mirroring of points X at A to yield Y
disp('check mirroring of points X at A to yield Y ...')
A = A/norm(A(1:3));
Normal = A(1:3);
S = -A(4);
H = [eye(3)-2*(Normal*Normal'), 2*S*Normal;[0,0,0,1]]
Xs = (H*[X,ones(3,1)]')';
disp('difference should be 0:')
check_mirroring = Xs(:,1:3)-Y

% projection
disp('projection ...')
disp('rotation matrix') % freely chosen
R = calc_Rot_r([0.1,0.2,-0.3])
disp('projection matrix') % translation freely chosen
P = R'*[eye(3), [2,2,-1.0]']

% image points
disp('projected image points')
x = (P*[X,ones(N,1)]')'
y = (P*[Y,ones(N,1)]')'

figure('Color','w')
plot(x(:,1)/x(:,3),x(:,2)/x(:,3),'or')
hold on
plot(y(:,1)/y(:,3),y(:,2)/y(:,3),'xb')
axis equal
title('Original points (red), mirrowed points (blue)')

disp(' ')
disp('Given image points, solve for E-matrix ... ')
%% solve for rotations
disp('solve for rotations ...')
% normal
n = cross(cross(x(2,:),y(2,:)),cross(x(3,:),y(3,:)))';
disp('normal')
n = n/norm(n)
disp('|[x1 cross y1, x2 cross y2, x3 cross y3]|, should be 0')
det_p = det([cross(x(1,:),y(1,:))',cross(x(2,:),y(2,:))',cross(x(3,:),y(3,:))'])

% Parameters
%% I. Solution
disp('1. Solution ....')
s = 1/(n(2)^2+n(3)^2)*(1-n(1))
q2 = -n(3)*s
q3 = n(2)*s

% Rotations
disp('rotations')
Rl = calc_Rot_q([1, 0,  q2,  q3])
Rr = calc_Rot_q([1, 0, -q2, -q3])
disp('Essential matrix')
E  = Rl*calc_S([1,0,0]')*Rr'

% check coplanarity constraints of point pairs
disp('check coplanarity constraints of point pairs ...')
disp('contradictions w_i = x_i^T E x_i')
w1 = x(1,:)*E*[-y(1,1),y(1,2),y(1,3)]'                                     %#ok<*NASGU>
w2 = x(2,:)*E*[-y(2,1),y(2,2),y(2,3)]'
w3 = x(3,:)*E*[-y(3,1),y(3,2),y(3,3)]'

% check whether points are in front of the cameras
% see lines 11-13 in Alg. 20, p. 583
N = cross(cross([1,0,0]',Rl'*x(1,:)'),[1,0,0]');
M = calc_S(N)* cross((Rl'*x(1,:)'),Rr'*[-y(1,1),y(1,2),y(1,3)]');
ssa = sign([1,0,0] * M);
sra = ssa.* sign((Rr'*[-y(1,1),y(1,2),y(1,3)]')'* N);
% if positive: in front
disp('positive, if points in front of camera (here: correct solution)')
[sra,ssa]

%% II. solution
disp(' ')
disp('2. Solution ....')
s = 1/(n(2)^2+n(3)^2)*(-1-n(1))
q2 = -n(3)*s
q3 = n(2)*s

% Rotations
disp('rotations')
Rl2 = calc_Rot_q([1, 0,  q2,  q3])
Rr2 = calc_Rot_q([1, 0, -q2, -q3])
disp('Essential matrix')
E  = Rl2*calc_S([1,0,0]')*Rr2'

% check coplanarity constraints of point pairs
disp('check coplanarity constraints of point pairs ...')
disp('contardictions w_i = x_i^T E x_i')
w1 = x(1,:)*E*[-y(1,1),y(1,2),y(1,3)]'
w2 = x(2,:)*E*[-y(2,1),y(2,2),y(2,3)]'
w3 = x(3,:)*E*[-y(3,1),y(3,2),y(3,3)]'

% check whether points are in front of the cameras
% see lines 11-13 in Alg. 20, p. 583
N = cross(cross([1,0,0]',Rl2'*x(1,:)'),[1,0,0]');
M = calc_S(N)* cross((Rl2'*x(1,:)'),Rr2'*[-y(1,1),y(1,2),y(1,3)]');
ssa = sign([1,0,0] * M);
sra = ssa.* sign((Rr2'*[-y(1,1),y(1,2),y(1,3)]')'* N);
% if positive: in front
disp('positive, if points in front of camera (here: wrong solution)')
[sra,ssa]

Rl'*Rl2
