%% test sugr_E_matrix
%
%
% Wolfgang Förstner 03/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

clearvars
close all

addpath(genpath('../../../General-Functions'));

%% 
disp('------------------')
disp('E = sugr_E_Matrix(A)')
A = randn(3);
E = sugr_E_Matrix(A);
normE = norm(E.E(:))                                                       %#ok<*NOPTS,*NASGU>
rankE = rank(E.E)


%%
disp('--------------------')
disp('E = sugr_E_Matrix(b,R)')
b = randn(3,1);
R = calc_Rot_r(randn(3,1));
E = sugr_E_Matrix(b,R);
bR = E.bR
normE = norm(E.E(:))
rankE = rank(E.E)
%% 
%
% 
disp('--------------------')
disp('E = sugr_E_matrix(b,R,CbRbR)')
b = randn(3,1);
R = calc_Rot_r(randn(3,1));
B = randn(5);
CbRbR = B'*B;
E = sugr_E_Matrix(b,R,CbRbR);
bR = E.bR
checkR = R'*R-eye(3)
normE = norm(E.E(:))
rankE = rank(E.E)
rankCee = rank(E.Cee)
checkCbRbR = E.CbRbR-pinv(E.JebR)*E.Cee*pinv(E.JebR)'

%
disp('--------------------')
disp('E=sugr_E_matrix(b,R,Cee)')
b = randn(3,1);
R = calc_Rot_r(randn(3,1));
B = randn(9);
Cee = B'*B;
E = sugr_E_Matrix(b,R,Cee);
bR = E.bR
checkR = R'*R-eye(3)
normE = norm(E.E(:))
rankE = rank(E.E)
rankCee = rank(E.Cee)
checkCee = E.Cee-E.JebR*E.CbRbR*E.JebR'



