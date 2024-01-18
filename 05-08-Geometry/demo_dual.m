%% DEMO Duality
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 07/17
% wenzel@igg.uni-bonn.de

clc
close all

addpath(genpath('../General-Functions'))

fprintf('\n------ DEMO Duality ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')
fprintf('\nGiven: 3D points\n')
X = random('bino',20,.5,[4,1]); disp(['X = [', num2str(X'),']'])
Y = random('bino',20,.5,[4,1]); disp(['Y = [', num2str(Y'),']'])
% any two auxiliary points  -> on planes through X and Y
U = random('bino',5,.5,[4,1]); disp(['U = [', num2str(U'),']'])
V = random('bino',5,.5,[4,1]); disp(['V = [', num2str(V'),']'])

%% plane through X and Y

fprintf('\nPlanes through X and Y\n')
A = null([X,Y,U]');
disp(['A = [', num2str(A'),']'])

A = calc_Gammadual(calc_Pi(X)*Y)'*U;
disp(['A = [', num2str(A'),']'])

B = null([X,Y,V]');
disp(['B = [', num2str(B'),']'])

B = calc_Gammadual(calc_Pi(X)*Y)'*V;
disp(['B = [', num2str(B'),']'])

%% check incidence
fprintf('\nCheck incidence:\n')
Incidence = [A'*X A'*Y; B'*X B'*Y];
disp(Incidence)

%% Lines
fprintf('Lines\n')
LXY = calc_Pi(X)*Y; disp(['L_XY = X cap Y = [', num2str(LXY'),']'])
LAB = calc_Pidual(A)*B; disp(['L_AB = A ^ B = [', num2str(LAB'),']'])

fprintf('normalize L_AB\n')
LAB = LAB/LAB(1)*LXY(1); disp(['L_AB = [', num2str(LAB'),']'])

%% check identity
fprintf('\nCheck identity\n')
should_be_zero = LXY'*calc_Dual*LAB;
disp(['Check should be zero:  ',  num2str(should_be_zero)])


