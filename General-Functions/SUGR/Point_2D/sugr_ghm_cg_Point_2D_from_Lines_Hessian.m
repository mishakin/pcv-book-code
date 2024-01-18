%% Determines residuals etc. for 2D point from passing lines 
% using line's Hessian representation
% 
% [va,C,cg,atr,btr] = sugr_ghm_cg_2d_point_from_lines_Hessian(l,le,xe,Cee)
%
% with Hessian representation
%
% * l   = 2 x 1 vector, observed line
% * le  = 2 x 1 vector, approximated fitted line
% * xe  = 2 x 1 vector, approximated point
% * Cee = 2 x 2 matrix, covariance matrix of observations
%
% * va  =   = 2 x 1 vector, approximate residual
% * cg  = scalar, residual of constraint
% * atr  = 1 x 2 vector, transposed Jacobian for xe
% * btr  = 1 x 2 vector, transposed jacobian for le
%
% Wolfgang Förstner 07/2012, adapted for sugr 6/2010
% wfoerstn@uni-bonn.de
% 
% See also sugr_Point_2D, sugr_Line_2D,  sugr_estimation_ml_Point_2D_from_Lines_Hessian

function [va,C,cg,atr,btr] = sugr_ghm_cg_Point_2D_from_Lines_Hessian(l,le,xe,C)

% g(l,x): x1 cos(l1) + x2 sin(l1) - l2 = 0

% approximate residual
va = le-l;
% Jacobian for x: 
atr = [cos(le(1)),sin(le(1))]; % row 2-vector

% Jacobain for l
btr = [-xe(1)*sin(le(1))+xe(2)*cos(le(1)),-1]; % row 2-vector
 
% residual of constraint
%cg = - l' * xe ; % scalar
cg = -(xe'*[cos(l(1));sin(l(1))]-l(2)); 
