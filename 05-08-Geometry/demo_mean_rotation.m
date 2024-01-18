%% DEMO  3D Mean Rotation
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 07/17
% wenzel@igg.uni-bonn.de

clc
close all

addpath(genpath('../General-Functions'))

fprintf('\n------ DEMO  3D Mean Rotation ------\n')
fprintf('\nIn order to shorten the output, we display transposed vectors\n')

% control of random numbers
sigma = 0.1;

% generate data
q_true = rand(4,1);
% Given: true 3D rotation
q_true = q_true/norm(q_true);
R_true = calc_Rot_q(q_true);

fprintf('\nGiven: noisy 3D rotations\n')
fprintf('\nQuaternions\n')
q1 = q_true + rand(4,1) * sigma; disp(['q1 = [', num2str(q1'),']'])
q2 = q_true*sign(rand) + rand(4,1) * sigma ; disp(['q2 = [', num2str(q2'),']'])
fprintf('\nAccording rotation matrices\n')
R1 = calc_Rot_q(q1)                                                        %#ok<NOPTS>
R2 = calc_Rot_q(q2)                                                        %#ok<NOPTS>

%% Average quaternions
fprintf('\nAverage quaternion\n')
M = q1 * q1' + q1 * q1';
[Ev,lambdas] = eig(M);
q_est = Ev(:,4); disp(['q_est = [', num2str(q_est'),']'])

fprintf('\nResidual\n')
dq = calc_Mq(q_true)\q_est;
difference_q_qtrue = norm(dq(2:4));
disp(['diff_q = ', num2str(difference_q_qtrue)])

%% Average rotation matrices
fprintf('\nAverage rotation matrices\n')
R_sum = R1 + R2;
[U,D,V] = svd(R_sum);
R_est = U * V'                                                             %#ok<NOPTS>

fprintf('\nResiduals\n')
difference_R_Rtrue = norm((R_est*R_true'-eye(3))/2);
disp(['diff_R = ', num2str(difference_R_Rtrue)])

fprintf('\nComparision: get rotation matrix from estimated quaternion\n')
R_true_from_q = calc_Rot_q(q_est)                                          %#ok<NOPTS>