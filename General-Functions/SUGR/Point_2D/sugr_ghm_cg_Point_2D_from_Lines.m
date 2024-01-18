%% Determines residuals etc. for 2D point from passing lines 
% 
% [lr,Crrt,cg,atr,btr] = sugr_ghm_cg_Point_2D_from_Lines(l,le,xe,Crro)
%
% in spherical mode with minimal representation
%
% * l   = 3 x 1 vector, observed line
% * le  = 3 x 1 vector, approximated fitted line
% * xe  = 3 x 1 vector, approximated point
% * Crro = 2 x 2 matrix, reduced covariance matrix of observations
%
% * lr  = 2 x 1 vector, reduced observation
% * Crrt = 2 x 2 matrix, reduced covariance matrix of transformed observations
% * cg  = scalar, residual of constraint
% * atr  = 1 x 2 vector, transposed Jacobian for xe
% * btr  = 1 x 2 vector, transposed jacobian for le
%
% Wolfgang Förstner 09/2010, adapted for sugr 1/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_Line_2D,  sugr_estimation_ml_Point_2D_from_Lines


function [lr,Cr,cg,atr,btr] = sugr_ghm_cg_Point_2D_from_Lines(l,le,xe,C)

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
atr = le' *  null(xe'); % row 2-vector

% Jacobain for l
btr = xe' * Jle; % row 2-vector
 
% residual of constraint
%cg = - l' * xe ; % scalar
cg = - le' * xe - btr * lr; 

