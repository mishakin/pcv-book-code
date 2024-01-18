%% exercise 12.9 conditioning
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

disp('-------------------------------------------------------------')
disp('Result of exercise 12.9 on conditioning the projection matrix')
disp('-------------------------------------------------------------')

% given projection matrix
P = [0.0012682807017 -0.0006478859649   0.0003109824561  -0.1611793859649
     0.0008165263157  0.0010670263157   0.0001048421052  -2.1127144736842 
     0.0000002017543  0.0000000350877   0.0000004561403  -0.0008359649122];

% transformation matrix for conditioning image coordinates
disp('Transformation matrix for conditioning image coordinates')
Tx = [1 0 -1500;0 1 -1000;0 0 1500]                                        %#ok<*NOPTS>

% transformation matrix for conditioning object coordinates
disp('Transformation matrix for conditioning object coordinates')
TX = [1 0 0 -500; 0 1 0 -1500;0 0 1 -150; 0 0 0 100]

% conditioned projection matrix
disp('Conditioned projection matrix')
Pb = (Tx)*P*inv(TX)                                                        %#ok<MINV>

% condition number of given projection matrix
disp('Condition number of given projection matrix')
condition_P = cond(P)

% condition number of conditioned projection matrix
disp('Condition number of conditioned projection matrix')
condition_Pb = cond(Pb)

disp('The routine ''cond'' determines the condition number')
disp('of a rectangular matrix from its singular values:')
singular_values = svd(Pb)
condition_number = max(singular_values)/min(singular_values)
