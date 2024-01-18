%% Demo 3D-geometry
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


fprintf('\n------ DEMO 3D projective geometry ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')
fprintf('\nGiven: 3D points\n')
X1 = [1,1,2,1]'; disp(['X_1 = [', num2str(X1'),']']);
X2 = [0,2,3,1]'; disp(['X_2 = [', num2str(X2'),']']);
X3 = [0,1,3,1]'; disp(['X_3 = [', num2str(X3'),']']);
X4 = [1,2,3,1]'; disp(['X_4 = [', num2str(X4'),']']);

fprintf('\nLine L_2 = X_3 ^ X_1 \n')
L2 = calc_Pi(X3)*X1; disp(['L_2 = [', num2str(L2'),']'])

fprintf('\nPlane A_4 = L_2 ^ X_2 \n')
A4 = calc_Gammadual(L2)*X2; disp(['A_4 = [', num2str(A4'),']'])

fprintf('\nGiven: Planes\n')
fprintf('A_1 parallel to XY-plane\n')
A1 = [0,0,-1,3]'; disp(['A_1 = [', num2str(A1'),']'])
fprintf('\nA_z = XY-plane\n')
Az=[0,0,1,0]'; disp(['A_z = [', num2str(Az'),']'])

fprintf('\nIntersection L_iz = A_1 cap A_z\n')
fprintf('\nLine at infinity perpendicular to Z-axis\n')
Linf_z = calc_Pidual(A1)*Az; disp(['Linf_z = [', num2str(Linf_z'),']'])