%% demo: check of signs single and two view geometry
%
% 1. Epipoles as projections and vectors of determinants 
% 2. Signs of epipolar lines 
% 3. Positivity of projection centre Z_h > 0
% 4. l' e'+l'' e'' = 0, PCV (12.286)
% 5. Sign of reconstructed 3D lines, PCV (13.288)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

clc

addpath(genpath('../General-Functions'))

disp('----------------------------------------------------------------')
disp('----- check Werner/Pajdla_s ls*es+lss*ess = 0 PCV (12.286) -----')
disp('----------------------------------------------------------------')

disp('Create Fundamental matrix F, assuming randomly given positive projection centres Z1 and Z2 ')
% random Projection matrices with positive projection centres
P1  = randn(3,4);                    % random P
DP1 = det(P1(1:3,1:3));              % determinant of A
P1  = P1/sign(DP1);                  % proper P, see PCV p. 474
                                     % --> positive projection centre                                  
% Projection centre as intersection of principal planes, see PCV (12.45)                                  
Z1  = -calc_Gammadual(calc_Pi(P1(1,:)')*P1(2,:)')'*P1(3,:)';

% same for P2
P2  = randn(3,4);
DP2 = det(P2(1:3,1:3));
P2  = P2/sign(DP2);
Z2  = -calc_Gammadual(calc_Pi(P2(1,:)')*P2(2,:)')'*P2(3,:)';

% corresponding Q's
Q1 = calc_Q_from_P(P1);             
Q2 = calc_Q_from_P(P2);

% Fundamental matrix, see PCV (13.19)
F = Q1*calc_Dual*Q2'                                                       %#ok<NOPTS>

% epipoles, see PCV (13.71)
e1 = P1*Z2;
e2 = P2*Z1;

% epipoles, see PCV (13.72)
e1_det = [...
    det([P2(1,:)',P2(2,:)',P2(3,:)',P1(1,:)']);...
    det([P2(1,:)',P2(2,:)',P2(3,:)',P1(2,:)']);...
    det([P2(1,:)',P2(2,:)',P2(3,:)',P1(3,:)']);...
      ];

disp(' ')
disp('----- Check epipoles -----')  
check_e1_minus_e1det_zero = (e1 - e1_det)';
disp('Epipole e1 = P1 times Z2 (13.71) minus e1 from camera planes in the ')
disp(['projection matrices (13.72) should be zero vector:   [', num2str(check_e1_minus_e1det_zero),']'])

check_e1F_zero=(e1'*F)';
check_Fe2_zero=(F*e2)';
disp(' ')
disp(['e1^T times F should be zero vector:   [', num2str(check_e1F_zero'),']^T'])
disp(['F times e2 should be zero vector:     [', num2str(check_Fe2_zero),']'])

disp(' ')
disp('----- Create random 3D Line and according projections into images ...')
% random 3D Line 
X1 = randn(4,1);            % two random 3D points
X2 = randn(4,1);
L  = calc_Pi(X1)*X2;        % resulting 3D line
l1 = Q1*L;                  % lines in images
l2 = Q2*L;

disp('Object point X1 is projected into image 1 as x11, into image 2 as x12')
disp('Object point X2 is projected into image 2 as x21, into image 2 as x22')
% projections
x11 = P1*X1;                % first point in images
x12 = P2*X1;
x21 = P1*X2;                % second point in images
x22 = P2*X2;

disp(' ')
disp('----- Calculate epipolar lines l1 and l2 of X1 ...')
% epipolar lines l1 and l2 of X1
l1x = cross(e1,x11)                                                        %#ok<NOPTS>
l1xF = F*x12;
l2x = cross(e2,x12)                                                        %#ok<NOPTS>
l2xF = F'*x11;

disp(' ')
disp('----- Check epipolar lines')
check_l1x_minus_l1xF_zero = (l1x-l1xF)';
check_l2x_minus_l2xF_zero = (l2x-l2xF)';
disp(['epipolar line l1 = e1 cross x11 minus l1 = F times x12 should be zero vector:     [', ...
    num2str(check_l1x_minus_l1xF_zero),']'])
disp(['epipolar line l2 = e2 cross x12 minus l2 = F^T times x11 should be zero vector:   [', ...
    num2str(check_l2x_minus_l2xF_zero),']'])

check_l_on_x_1_zero = (l1x'*x11)';
disp(['image point x11 should be on l1, thus l1^T times x11 should be zero:              ', num2str(check_l_on_x_1_zero)])
check_e_on_x_1_zero = (l1x'*e1)';
disp(['epipole e1 should be on l1, thus l1^T times e1 should be zero:                    ', num2str(check_e_on_x_1_zero)])

check_l_on_x_2_zero = (l2x'*x12)';
disp(['image point x12 should be on l2, thus l2^T times x12 should be zero:              ', num2str(check_l_on_x_2_zero)])
check_e_on_x_2_zero = (l2x'*e2)';
disp(['epipole e2 should be on l2, thus l2^T times e2 should be zero:                    ', num2str(check_e_on_x_2_zero)])

%% start of check
disp(' ')
disp('---- Check positive sign and value of Z0 vs abs|det(A)| , see PCV p. 474 -----')
Z_s_positive_Dets       = [Z1',abs(DP1);Z2',abs(DP2)]                      %#ok<NOPTS>


disp('----- Check sign of l1*e1 and l2*e2 -----')
lses   = l1'*e1;
lssess = l2'*e2;
check_lses_plus_lssess_zero = lses+lssess;
disp(['PCV (13.287)   l1*e1 + l2*e2 should be zero:        ', num2str(check_lses_plus_lssess_zero)])
check_QP_minus_PiZ_zero = calc_Dual*Q1'*P1-calc_Pi(Z1);
disp(['PCV (13.288)   Q_dual(P1) - Pi(Z1) should be zero:  ', num2str(check_lses_plus_lssess_zero)])
disp(' ')
%% 
disp('Check sign of 3D line reconstruction PCV (13.288)')
A1 = P1'*l1;
A2 = P2'*l2;
Lc = sign(l1'*e1)*calc_Pidual(A1)*A2;   
disp( ['L = A2 cap A1 =                                    [',num2str(Lc'),']'])
% Lc/Lc(1)*sign(Lc(1));
% L/L(1)*sign(L(1));
factor_12 = sign(l1'*e1);
disp(['Sign of 3D line: ', num2str(factor_12)])
check_L_and_Lreconstructed_zero = (Lc/Lc(1)*sign(Lc(1))-L/L(1)*sign(L(1)))';
disp(['Diff L and reconstructed 3D line should be zero:   [', num2str(check_L_and_Lreconstructed_zero),']'])

disp(' ')
Lc = sign(l2'*e2)*calc_Pidual(A2)*A1;
disp( ['L = A1 cap A2 =                                    [',num2str(Lc'),']'])
% Lc/Lc(1)*sign(Lc(1));
% L/L(1)*sign(L(1));
factor_21 = sign(l2'*e2);
disp(['Sign of 3D line: ', num2str(factor_21)])
check_zero = (Lc/Lc(1)*sign(Lc(1))-L/L(1)*sign(L(1)))';
disp(['Diff L and reconstructed 3D line should be zero:   [', num2str(check_zero),']'])





