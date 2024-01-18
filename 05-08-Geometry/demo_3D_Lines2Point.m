%% Demo intersection and join of 3D lines
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 07/17
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

addpath(genpath('../General-Functions'))

fprintf('\n------ Demo point and plane closest to two 3D lines ------\n')

% Given points on lines
X = [ 2,1,2,1]';
Y = [ 6,1,2,1]';
% intersecting lines
Z = [ 3,2,2,1]';
T = [ 3,7,2,1]';
% % non-interesecting lines
% Z = [ 3,2,6,1]';
% T = [ 3,7,6,1]';
% line at infinity
% Z = [ 3,2,2,0]';
% T = [ 3,7,2,0]';

fprintf('\nIn order to shorten the output, we display transposed vectors\n')


fprintf('\nGiven: 3D points\n')

disp(['X = [', num2str(X'),']'])
disp(['Y = [', num2str(Y'),']'])
disp(['Z = [', num2str(Z'),']'])
disp(['T = [', num2str(T'),']'])
fprintf('\nGiven: 3D lines\n')
L = calc_Pi(X)*Y; disp(['L = X ^ Y = [', num2str(L'),']'])

M = calc_Pi(Z)*T; disp(['M = Z ^ T = [', num2str(M'),']'])


%% Intersection point X of two 3D-lines L and M
fprintf('\nPoint closest to lines\n')

% general solution: any position but lines are not equal
Intersection_general = calc_intersect_3D_Lines2Point(L,M);
disp(    ['algebraic: X(L,M)             = [', num2str(Intersection_general'),']'])
if Intersection_general(4) ~= 0
    disp(['algebraic: X(L,M), normalized = [', num2str(Intersection_general'/Intersection_general(4)),']'])
end

% finite solution: for lines that are not parallel
Intersection_finite = calc_intersect_3D_Lines2Point_finite(L,M);
disp(['finite: X(L,M)                = [', num2str(Intersection_finite'),']'])
if Intersection_finite(4) ~= 0
    disp(['finite: X(L,M), normalized    = [', num2str(Intersection_finite'/Intersection_finite(4)),']'])
end

%% Plane through two 3D-lines L and M
fprintf('\nPlane closest to lines\n')

% general solution: lines are coplanar and not equal
Plane_general = calc_join_3D_Lines2Plane(L,M);
disp(['algebraic: A(L,M)             = [', num2str(Plane_general'),']'])
if norm(Plane_general(1:3)) ~= 0
    disp(['algebraic: A(L,M), normalized = [', num2str(Plane_general'/norm(Plane_general(1:3))),']'])
end

% finite solution: for lines that are not parallel
Plane_finite = calc_join_3D_Lines2Plane_finite(L,M);
disp(['finite: A(L,M)                = [', num2str(Plane_finite'),']'])
if norm(Plane_finite(1:3)) ~= 0
    disp(['finite: A(L,M), normalized    = [', num2str(Plane_finite'/norm(Plane_finite(1:3))),']'])
end