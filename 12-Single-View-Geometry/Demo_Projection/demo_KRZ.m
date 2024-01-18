%% Demo partitioning of P into K, R, and Z
% 
% PCV, Algorithm 18, p. 500 
%
% The result should not change if the sign of P changes.
% The sign of the focal length c should be the required one
%
% vary:
% - sign of K(3,3), changing sign of P
% - sign of c
%
% check the correctnes of the decomposition for all four cases
% diag(K) = [+1,+1,+1] : require c = +1
% diag(K) = [-1,-1,+1] : require c = -1
% diag(K) = [+1,+1,-1] : require c = -1
% diag(K) = [-1,-1,-1] : require c = +1
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 04/18
% wenzel@igg.uni-bonn.de

clc
close all
clearvars

addpath(genpath('../../General-Functions'))

fprintf('\n------ Demo partitioning of P into K, R, and Z ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')

disp(' ')
% R = eye(3);
disp('Common rotation matrix')
R = calc_Rot_q([1,.0,0.0,2])                                               %#ok<NOPTS>

%% 1. camera
fprintf('\n-------- camera 1: K(3,3) > 0, c > 0  --------\n')
disp('create P with random Z and')
K = [+1,0.01,0.01;0,+1,0.1;0,0,1];
disp(['diag(K) = [', num2str(diag(K)'),']'])

Z = rand(3,1); %zeros(3,1);
P = calc_P_from_KRZ(K,R,Z)                                                 %#ok<NOPTS>
d = calc_viewing_direction(P);
disp(['Viewing direction d = [', num2str(d'),']'])

disp(' ')
disp('Get K, R, Z from P with positive c')
[K_,R_,Z_,~] = calc_KRZ_from_P(P,1);
K_                                                                         %#ok<NOPTS>
R_                                                                         %#ok<NOPTS>
Z_'                                                                        %#ok<NOPTS>

disp('Check calculated values')
dK1 = K/K(3,3)-K_;
disp(['norm(dK) = [', num2str(norm(dK1)),']'])
dR1 = R-R_;
disp(['norm(dR) = [', num2str(norm(dR1)),']'])
dZ1 = Z-Z_;
disp(['norm(dZ) = [', num2str(norm(dZ1)),']'])

%% 2. camera
fprintf('\n-------- camera 2: K(3,3) > 0, c < 0 --------\n')
disp('Again, create P with random Z and')
K=[-1,0.01,0.01;0,-1,0.1;0,0,1];
disp(['diag(K) = [', num2str(diag(K)'),']'])

Z = rand(3,1); % zeros(3,1);
P = calc_P_from_KRZ(K,R,Z)                                                 %#ok<NOPTS>
d = calc_viewing_direction(P);
disp(['Viewing direction d = [', num2str(d'),']'])

disp(' ')
disp('Get K, R, Z from P with negative c')
[K_,R_,Z_,~] = calc_KRZ_from_P(P,-1);
K_                                                                         %#ok<NOPTS>
R_                                                                         %#ok<NOPTS>
Z_'                                                                        %#ok<NOPTS>

disp('Check calculated values')
dK2 = K/K(3,3)-K_;
disp(['norm(dK) = [', num2str(norm(dK2)),']'])
dR2 = R-R_;
disp(['norm(dR) = [', num2str(norm(dR2)),']'])
dZ2 = Z-Z_;
disp(['norm(dZ) = [', num2str(norm(dZ2)),']'])


%% 3. camera
fprintf('\n-------- camera 3: K(3,3) < 0, c < 0 --------\n')
disp('Again, create P with random Z and')
K=[+1,0.01,0.01;0,+1,0.1;0,0,-1];
disp(['diag(K) = [', num2str(diag(K)'),']'])

Z = rand(3,1); % zeros(3,1);
P = calc_P_from_KRZ(K,R,Z)                                                 %#ok<NOPTS>
d = calc_viewing_direction(P);
disp(['Viewing direction d = [', num2str(d'),']'])

disp(' ')
disp('Get K, R, Z from P with negative c')
[K_,R_,Z_,~] = calc_KRZ_from_P(P,-1);
K_                                                                         %#ok<NOPTS>
R_                                                                         %#ok<NOPTS>
Z_'                                                                        %#ok<NOPTS>

disp('Check calculated values')
dK3 = K/K(3,3)-K_;
disp(['norm(dK) = [', num2str(norm(dK3)),']'])
dR3 = R-R_;
disp(['norm(dR) = [', num2str(norm(dR3)),']'])
dZ3 = Z-Z_;
disp(['norm(dZ) = [', num2str(norm(dZ3)),']'])

%% 4. camera
fprintf('\n-------- camera 4: K(3,3) < 0, c > 0 --------\n')
disp('Again, create P with random Z and')
K=[-1,0.01,0.01;0,-1,0.1;0,0,-1];
disp(['diag(K) = [', num2str(diag(K)'),']'])

Z = rand(3,1); % zeros(3,1);
P = calc_P_from_KRZ(K,R,Z)                                                 %#ok<NOPTS>
d = calc_viewing_direction(P);
disp(['Viewing direction d = [', num2str(d'),']'])

disp(' ')
disp('Get K, R, Z from P with positive c')
[K_,R_,Z_,~] = calc_KRZ_from_P(P,1);
K_                                                                         %#ok<NOPTS>
R_                                                                         %#ok<NOPTS>
Z_'                                                                        %#ok<NOPTS>

disp('Check calculated values')
dK4 = K/K(3,3)-K_;
disp(['norm(dK) = [', num2str(norm(dK4)),']'])
dR4 = R-R_;
disp(['norm(dR) = [', num2str(norm(dR4)),']'])
dZ4 = Z-Z_;
disp(['norm(dZ) = [', num2str(norm(dZ4)),']'])

sum_check=...
    norm(dK1(:))+norm(dK2(:))+norm(dK3(:))+norm(dK4(:))+...
    norm(dR1(:))+norm(dR2(:))+norm(dR3(:))+norm(dR4(:))+...
    norm(dZ1(:))+norm(dZ2(:))+norm(dZ3(:))+norm(dZ4(:));

disp(' ')
display(['Total sum of absolute differences =',num2str(sum_check)]);
