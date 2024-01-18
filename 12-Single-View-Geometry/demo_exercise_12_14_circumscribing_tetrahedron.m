%% Demo Barycentric coordinates of point cloud within minimal spanning
% tetrahedron
%
% see Kaseorg, A. (2014). 
% How do I find the side of the largest cube completely contained inside a 
% regular % tetrahedron of side s? 
% https://www.quora.com/How-do-I-find..., last visited 10/2015. 538
%
% see PCV exercise 12.14
%
% Wolfgang Förstner 6/2016
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

disp('--------------------------------------------')
disp('   Barycentric coordinates of point cloud   ')
disp('--------------------------------------------')

% set range of coordinates
% ddx=5;ddy=20;ddz=0.1;
ddx=1;ddy=1;ddz=1;

cube=0; % Cube given by params dy, dy, dz or unitcube

if cube == 1
    disp('Cube [-1,+1]^3')
    X=[[-1,1,1,-1,-1,1,1,-1]*ddx;
       [1,1,-1,-1,1,1,-1,-1]*ddy;
       [-1,-1,-1,-1,1,1,1,1]*ddz];
else
    disp(['Cube ',num2str(ddx),' x ',num2str(ddy),' x ',num2str(ddz)])
    N = 1000;
    X = [rand(1,N)*ddx;rand(1,N)*ddy;rand(1,N)*ddz];
end

% determine minimal bounding box
xmin = min(X(1,:));
ymin = min(X(2,:));
zmin = min(X(3,:));
xmax = max(X(1,:));
ymax = max(X(2,:));
zmax = max(X(3,:));
mx = ( xmax + xmin )/2;
my = ( ymax + ymin )/2;
mz = ( zmax + zmin )/2;

% four points of minimal spanning tetrahedron
dx = xmax - xmin;
dy = ymax - ymin;
dz = zmax - zmin;
u = ( 2+sqrt(2) )/6;
v = u/sqrt(2);
T = [mx-u*dx, my, mz-v*dz;...
     mx+u*dx, my, mz-v*dz;...
     mx, my-u*dy, mz+v*dz;...
     mx, my+u*dy, mz+v*dz]';
disp('Corners of Tetrahedron')
disp(num2str(T));

Th = [T;ones(1,4)];

% check barycentric coordinates
alpha = inv(Th)*[X;ones(1,N)];                                             %#ok<MINV>

% max alpha should be < 1
disp(['maximum absolute Barycentric coordinate =',num2str(max((alpha(:))))])

figure
hold on
plot3(X(1,:),X(2,:),X(3,:),'.b')
plot3(  [ T(1,1), T(1,2), T(1,3), T(1,4), T(1,1), T(1,3), T(1,2), T(1,4)],...
        [ T(2,1), T(2,2), T(2,3), T(2,4), T(2,1), T(2,3), T(2,2), T(2,4)],...
        [ T(3,1), T(3,2), T(3,3), T(3,4), T(3,1), T(3,3), T(3,2), T(3,4)],...
        '-k','LineWidth',4' );

