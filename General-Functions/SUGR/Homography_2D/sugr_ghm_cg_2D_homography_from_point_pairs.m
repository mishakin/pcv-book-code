%% determines residuals etc. for 2D homography from point pairs
% with minimal representation
% 
% [lr,cg,Crrt,atr,btr] = sugr_ghm_cg_2D_homography_from_point_pairs(l,le,xe,Crro)
%
% * l   = 6 x 1 vector, observed point pair
% * le  = 6 x 1 vector, approximated fitted point pair
% * xe  = 3 x 3 matrix, approximated homography
% * Crro = 4 x 4 matrix, reduced covariance matrix of point pair
%
% * lr  = 4 x 1 vector, reduced observation
% * cg  = 2 x 1 vector, residual of constraint
% * Crrt = 4 x 4 matrix, reduced covariance matrix of transformed observations
% * atr  = 2 x 8 vector, transposed Jacobian for cg -> xe
% * btr  = 2 x 4 vector, transposed jacobian for cg -> le
%
% 0 != S^(l)(y)  H  x = -S^(l)(H x)  y = (x'  cron  S^(l)(y))  vec(H)
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 03/2011

function [lr,Cr,cg,atr,btr] = sugr_ghm_cg_2D_homography_from_point_pairs(l,le,xe,C)

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

% trasnferred covariance matrix
Cr = JR * C * JR';

% determine selected constraints
[Ssy, Rsy] = calc_S_reduced(le(4:6));

% Jacobian for cg -> x
Jhr = sugr_get_Jacobian_Jhr_Homography_2D(xe);
atr = kron(le(1:3)', Ssy) * Jhr; % 2 x 8 Matrix

% Jacobain for cg -> l
btr = [Ssy * xe * null(le(1:3)') , ...
      -Rsy * calc_S(xe * le(1:3)) * null(le(4:6)')] ; % 2 x 4 matrix
 
% residual of constraint (lr = -vr)
cg = - Ssy * xe * le(1:3)  - btr * lr; 

