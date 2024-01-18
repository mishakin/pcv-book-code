%% calculate line from two statistically independent points
%
% L = sugr_construct_join_line_3D(X,Y)
%
% X,Y    given 3D points as sugr objects
%
% L      joining 3D line as sugr object
% 
% Wolfgang Förstner 2013-09-17
% wfoerstn@uni-bonn.de 
%
% wf 2016-09-02 added variance propagation
% sw 9/2016
%
% See also sugr_Line_3D, sugr_minimal_3D_Line, sugr_constrain_3D_Line

function L = sugr_construct_join_Line_3D(X,Y)

%% determine homogeneous vectors

% Jacobians dL/dX, dL/dY
J_LX =   sugr_calc_Pi(X.h);
J_LY = - sugr_calc_Pi(Y.h);

h = J_LX*Y.h;

%% determine covariance matrix
% CovM of homogeneous 4-vectors
CovM_X = sugr_get_CovM_homogeneous_Vector(X);
CovM_Y = sugr_get_CovM_homogeneous_Vector(Y);

% covariance matrix of Plücker vector
CovM_Lh = J_LX * CovM_X *  J_LX' + J_LY * CovM_Y *  J_LY';

% generate sugr object 
L = sugr_Line_3D(h,CovM_Lh);
