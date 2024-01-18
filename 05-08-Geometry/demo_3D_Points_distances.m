%% demo distanc 3D point
%
% Susanne Wenzel 12/17
% wenzel@igg.uni-bonn.de

clc
close all
clearvars


addpath(genpath('../General-Functions'))

ss = plot_init;

fprintf('\n------ DEMO 3D point distances ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')

fprintf('\nGiven: 3D points\n')
X = [1,1,1]'; disp(['X = [', num2str(X'),']'])
X_ = [X;1];
Y = [3,1,1]'; disp(['Y = [', num2str(Y'),']'])
Y_ = [Y;1];
Z = [4,4,1]'; disp(['Z = [', num2str(Z'),']'])
Z_ = [Z;1];

figure('Color','w','Position',[50 ss(2)/3 ss(1)/3 ss(2)/2]);
pltp = @(X)scatter3(X(1),X(2),X(3),'ko');
plttxt = @(X,txt)text(X(1)+.1,X(2)-.1,X(3),txt);
hold on
pltp(X),plttxt(X,'X')
pltp(Y),plttxt(Y,'Y')
pltp(Z),plttxt(Z,'Z')
fill3([X(1),Y(1),Z(1)],[X(2),Y(2),Z(2)],[X(3),Y(3),Z(3)],'r')
xlim([0,5]),ylim([0,5]), zlim([0,5])
set(gca,'CameraPosition',[ -4.3227  -28.5303   31.9206])
set(gca,'CameraViewAngle',10.1830)

fprintf('\n3D line through X-Y\n')
L_ = calc_Pi(X_)*Y_; disp(['L = X ^ Y = [', num2str(L_'),']'])


fprintf('\nDistance of Z to line L\n')
d_ZL = calc_distance_3D_point_from_3D_line(Z_,L_)                           %#ok<*NOPTS>

fprintf('\nDistance of Z to linesegment X-Y\n')
d_Z_XY = calc_distance_3D_point_from_linesegment(Z,X,Y)

fprintf('\nDefine more 3D points\n\n')
T = [2.5  2.  4; ... % a point above the triangle
    2   0.5  0.5; ... % outside the triangle, next to one side
    0.5 0.5 0.5;... % outside the trinagle, next to one corner
    ]';

for i = 1:size(T,2)
    disp(['T_',num2str(i),' = [', num2str(T(:,i)'),']'])
    pltp(T(:,i)),plttxt(T(:,i),['$T_',num2str(i),'$'])
    
    d_T_XYZ = calc_distance_3D_point_from_triangle(T(:,i),X,Y,Z);
    disp(['Distance of T_',num2str(i),' to line L',...
        num2str(d_T_XYZ)])
    
    disp('press F5 to continue')
    keyboard
end