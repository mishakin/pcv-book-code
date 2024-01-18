%% determines residuals etc. for estimating plane from points
% 
% [lr,cg,Crrt,atr,btr] = sugr_ghm_cg_Plane_from_points(l,le,xe,Crro)
%
% in spherical mode with minimal representation
%
% * l   = 4 x 1 vector, observed point
% * le  = 4 x 1 vector, approximated fitted point
% * xe  = 4 x 1 vector, approximated plane
% * Crro = 3 x 3 matrix, reduced covariance matrix of observations
%
% * lr  = 3 x 1 vector, reduced observation
% * cg  = scalar, residual of constraint
% * Crrt = 3 x 3 matrix, reduced covariance matrix of transformed observations
% * atr  = 1 x 3 vector, transposed Jacobian for xe
% * btr  = 1 x 3 vector, transposed jacobian for le
%
% Wolfgang Förstner 07/2012, adapted for sugr 1/2011
% wfoerstn@uni-bonn.de 
%

function [lr,Cr,cg,atr,btr] = sugr_ghm_cg_Plane_from_points(l,le,xe,C)

% Jacobian for observation
Jle = null(le');

% reduced observation
lr  = Jle' * l;

% Rotation from l to le
R = calc_Rot_ab(le,l);

% Jacobian for Cr
JR = Jle' * R' * null(l');

% reduced covariance matrix
Cr = JR * C * JR';

% Jacobian for x
atr = le' *  null(xe'); % row 3-vector

% Jacobain for l
btr = xe' * Jle; % row 3-vector
 
% residual of constraint
%cg = - l' * xe ; % scalar
cg = - le' * xe - btr * lr; 

