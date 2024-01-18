%% Determines residuals etc. for projection from point pairs
% 
% [lr,cg,Crrt,atr,btr] = sugr_ghm_cg_Projection_3D_2D_from_point_pairs_e(l,le,xe,Crro)
%
% with euclidean representation
%
% * l   = 6 x 1 vector, observed point pair
% * le  = 6 x 1 vector, approximated fitted point pair
% * xe  = 3 x 4 matrix, approximated projection
% * Crro = 5 x 5 matrix, reduced covariance matrix of point pair
%
% * lr  = 5 x 1 vector, reduced observation
% * cg  = 2 x 1 vector, residual of constraint
% * Crrt = 5 x 5 matrix, reduced covariance matrix of transformed observations
% * atr  = 2 x 11 matrix, transposed Jacobian for cg -> xe
% * btr  = 2 x 5 matrix, transposed jacobian for cg -> le
%
% 0 != y - c(P  X) = y - c( (X' kron I3) vec(P)) ; c(x) = x(1:2)/x(3)
%
% Wolfgang Förstner 02/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs
% sugr_estimation_ml_Projection_3D_2D_from_point_pairs
% sugr_ghm_cg_Projection_3D_2D_from_point_pairs

function [lr,Cr,cg,atr,btr] = sugr_ghm_cg_Projection_3D_2D_from_point_pairs(l,le,xe,C)

%% preparation
% current estimate for projection matrix
Pe  = reshape(xe,3,4);

% homogeneous coordinates of projection
yh = Pe * le(1:4);


% 2x3 Jacobian of c(.) from homogenoues to Euclidean (12.129)
Jc = [yh(3)*eye(2), -yh(1:2)]/yh(3)^2; 

% Jacobian for observation -> reduced observation: 6 x 5 matrix
Jle = [null(le(1:4)')  zeros(4,2);... 
       zeros(2,3)      eye(2)] ;

% reduced observation_ 5 x 1 vector
lr  = [null(le(1:4)')' * l(1:4);  l(5:6)-le(5:6)];


%% adaptation of Covariance matrix
% Rotation from l to le: 6 x 6 matrix
R = [calc_Rot_ab(le(1:4),l(1:4)) zeros(4,2);...
    zeros(2,4) eye(2)];
% Jacobian for Cr
JR = Jle' * R' * Jle; %[null(le(1:4)')  zeros(4,2);... 
                      % zeros(2,3)      eye(2)] ;
% transferred covariance matrix
Cr = JR * C * JR';

%% Jacobians
% Jacobian for cg -> x
atr = - Jc *  kron(le(1:4)', eye(3)) * null(xe'); % 2 x 11 Matrix

% Jacobain for cg -> l
btr = [- Jc * Pe * null(le(1:4)') ,  eye(2)];  % 2 x 5 matrix
 
%% residual of constraint (lr = -vr)

% constraint = -(l(5:6)-yh(1:2)/yh(3));
cg = -  (le(5:6)- yh(1:2)/yh(3)) - btr * lr;
%disp(strcat('-cg(le,xe) = ',num2str(le(5:6)'- (yh(1:2)/yh(3))'),', -btr_lr =',num2str((btr * lr)')));

