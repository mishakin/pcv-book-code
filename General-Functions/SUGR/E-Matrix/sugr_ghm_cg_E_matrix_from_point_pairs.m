%% Constraints and Jacobians for GHM for estimating E-matrix from point pairs
% i.e., determines residuals etc. for E-matrix from point pairs
% with minimal representation
% 
% [lr,cg,Crrt,atr,btr] = sugr_ghm_cg_E_matrix_from_point_pairs(l,le,xe,Crro)
%                        
%
% * l   = 6 x 1 vector, observed point pair
% * le  = 6 x 1 vector, approximated fitted point pair
% * xe  = 3 x 4 matrix [b^a, R^a] of approximated values
% * Crro = 4 x 4 matrix, reduced covariance matrix of point pair
%
% * lr  = 4 x 1 vector, reduced observation
% * cg  = 1 x 1 vector, residual of constraint
% * Crrt = 4 x 4 matrix, reduced covariance matrix of transformed observations
% * atr  = 1 x 5 vector, transposed Jacobian for cg -> xe
% * btr  = 1 x 4 vector, transposed jacobian for cg -> le
%
% 0 != x'^T S(b) R' x'' = x''^T R S^T(b) x' = 
%
% Wolfgang Förstner 03/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

function [lr,Cr,cg,atr,btr] = sugr_ghm_cg_E_matrix_from_point_pairs(l,le,xe,C)

% approximate values
ba = xe(:,1);
Ra = xe(:,2:4);
Ea = calc_S(ba)*Ra';

% Jacobian for observation -> reduced observation: 6 x 4 matrix
Jle = [null(le(1:3)')      zeros(3,2);... 
       zeros(3,2)      null(le(4:6)')] ;

% reduced observation_ 4 x 1 vector
lr  = Jle' * l;

% Rotation from l to le: 6 x 6 matrix
R = [calc_Rot_ab(le(1:3),l(1:3)) zeros(3);...
    zeros(3) calc_Rot_ab(le(4:6),l(4:6))];

% Jacobian for Cr
JR = Jle' * R' * [null(le(1:3)')      zeros(3,2);... 
                  zeros(3,2)      null(le(4:6)')] ;

% transferred covariance matrix
Cr = JR * C * JR';


% Jacobian for cg -> x
xais    = le(1:3);
xaiss   = le(4:6);
xaiss1  = Ra' * xaiss;
lais    = Ea  * xaiss;
laiss   = Ea' * xais;
atr = [cross(xaiss1 , xais)' * null(ba') , cross(laiss , xaiss)']; % 1 x 5 vector

% Jacobain for cg -> l
btr = [lais' * null(xais') , laiss' * null(xaiss')] ; % 1 x 4 vector
 
% residual of constraint (lr = -vr)
cg = - le(1:3)'*calc_S(ba)*Ra'*le(4:6)  - btr * lr; 

