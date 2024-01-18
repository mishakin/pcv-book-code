%% demo: bR from (U,V,xs,xss)
%
% Check correctnes of basis and Rotation from 
% (U,V) or E and corresponding points (see Alg. 20)
%
% U,V      = from svd(E)
% xs,xss   = points in left and right camera system
% 
% bR       = [b,R] matrix of relative orientation 
%
% Wolfgang Förstner 6/2017
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

clc
clear vars

addpath(genpath('../../General-Functions'))

sugr_INIT

disp('--------------------------------------------')
disp(' bR from (U,V) = svd(E),x'',x'''' (Alg. 20) ')
disp('--------------------------------------------')

%% generate random E-matrix
% choose rotation matrix
Rs  = eye(3);
Rss = calc_Rot_r(randn(3,1)); %diag([1,-1,-1]);
% choose basis direction
bh  = randn(3,1); 
bh  = bh/norm(bh); %[-1,0,0]';
% choose 3d point
X = randn(3,1);
% determine image coordinates
xs  = Rs*X;
xss = Rss*(X-bh);

% choose E-matrix with random sign
v = sign(randn(1));
disp(strcat('random_sign of_E =  ',num2str(v)));
E = Rs*calc_S(bh)*Rss'*v                                                     %#ok<*NOPTS>

% provide SVD with positive U and V
[U,D,V] = svd(E);
U = U*det(U);
V = V*det(V);

% determine basis, Rotation and case for relative orientation
% (arguments xs' and xss' must be a list of row vectors, here only one)

%% from (E,xs,xss)
disp('b and R from E using set of pairs of directions')
% Alg. 20
[bes,Res] = sugr_bR_from_E_uv(E,xs',xss')

Es_check = sugr_E_Matrix(bes,Res);
disp('check results, differences of given and estimated b and R should be zero' )
diff_bRs_should_be_zero = [bh-bes,Rss-Res,v*E-Es_check.E]

%% from (U,V,xs,xss)
disp('b and R from U and V using set of pairs of directions')
% similar to Alg. 20, but starting from (U,V) instead from E
[be,Re] = sugr_select_bR_from_E_UV_uv(U,V,xs',xss')                              

E_check = sugr_E_Matrix(be,Re);
disp('check results, differences of given and estimated b and R should be zero' )
diff_bR_should_be_zero = [bh-be,Rss-Re,v*E-E_check.E]                        



