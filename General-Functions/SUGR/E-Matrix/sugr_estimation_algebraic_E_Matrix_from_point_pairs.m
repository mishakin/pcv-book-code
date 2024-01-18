%% Algebraic estimation of E-Matrix from point pairs, 
% assumption: independent observations
% assumption: points are conditioned.
%
% following Sect. 13.3.2.3, p. 575
%    with CovM following Sect. 4.9.2.4
% 
% [E,est_sigma_0,error] = sugr_estimation_algebraic_E_Matrix_from_point_pairs(P) 
%
% P contains point pairs
% P.h = [X,Y]
% * X = N x 3 matrix of points 
% * Y = N x 3 matrix of points
% P.Crr = N x 4 x 4 contains reduced covariance matrix
% 
% model
% X(n,:) * E * Y(n,:)' = 0
%
% E      relatiove orientation
%        E.bR = estimated homography with uncertainty, 
%               using algebraic minimization, thus neglecting the accuracy
%        E.Crr = CovM of reduced parameters, derived by variance propagation
%
% est_sigma_0  estimated standard deviation of constraint
%
% error   = 0  ok
%         = 1  no basis found
%
% Wolfgang Förstner 2/2010
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix, sugr_estimation_ml_E_Matrix_from_point_pairs


function [E,est_sigma_0,error] = sugr_estimation_algebraic_E_Matrix_from_point_pairs(P)

% pair coordinates
Ph   = P.h;
% homogeneous coordinates
xl = Ph(:,1:3);
xr = Ph(:,4:6);
% number of point pairs
N = size(Ph,1);

%% 
% estimate E algebraically using > 7 points
%
% build coefficient matrix
A = zeros(N,9);
for n = 1:N 
    % Kronecker 1x9 -matrix
    A(n,:)  = [xr(n,1)*xl(n,:), xr(n,2)*xl(n,:) xr(n,3)*xl(n,:)];
end

% partition
N_rows = size(A,1);
if N_rows == 8
    e = null(A);
    A_plus = A'*inv(A*A');                                                 %#ok<*MINV>
    est_sigma_0 = 1;
else
    [U,D,V] = svd(A,'econ');
    % find algebraically best solution
    e        = V(:,9);
    % determine A+ = pseudo inverse of A with last singular value = 0
    Di          = inv(D+eps*eye(9));
    est_sigma_0 = sqrt(Di(9,9)/(N-5));
    Di(9,9)  = 0;
    A_plus   = V * Di * U';
end
E0   = reshape(e,3,3);

% Enforce same eigenvalues == 1
[U0,~,V0]= svd(E0);
U0        = U0*det(U0);
V0        = V0*det(V0);
Ee        = U0 * diag([1,1,0]) * V0';

% select from four cases
[b,R] = sugr_select_bR_from_E_UV_uv(U0,V0,xl,xr);
if norm(b) == 0
    E           = zeros(3);
    est_sigma_0 = 0;
    error       = 1;
    return
end

% determine CovM of residuals of constraints from x=xl, y=xr
% (following Sect. 4.9.2.4)
% x' E y = g
% d(x' E y) = d(g)
% x' E dy + y' E' dx = dg
%x' E Jyr dyr + y' E' Jxr dxr = dg
%
Cgg = sparse(zeros(N));
for n=1:N
    x1    = xl(n,:)';
    x2    = xr(n,:)';
    % VarProp for g 
    J_x = + x2' * Ee' * null(x1');
    J_y = + x1' * Ee  * null(x2');
    J = [J_x, J_y];
    % effect of both
    Cgg(n,n)  = Cgg(n,n) ...
        + J * squeeze(P.Crr(n,:,:)) * J'; %#ok<SPRIX>
end

% determine Crr for 5 parameters
Chh = A_plus * Cgg * A_plus';

E     = sugr_E_Matrix(b,R,Chh);
error = 0;
