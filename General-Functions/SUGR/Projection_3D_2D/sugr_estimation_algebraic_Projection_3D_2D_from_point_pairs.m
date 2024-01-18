%% Algebraic estimate of Projection from corresponding 3D-2D point pairs
%
% assuming independent observations
% 
% P = sugr_estimation_algebraic_Projection_3D_2D_from_points(X,y) 
%
% X.e   = N x 3 matrix if points 
% y.e   = N x 2 matrix if points
% X.Cee = N x 3 x 3 matrix of covariances for 3D points
% y.Cee = N x 2 x 2 matrix of covariances for 2D points
% 
% model with constraints assuming points not at infinity
% 0 != Rm * S(y) * P * X                Rm = [eye(2), zeros(2,1)]
%    = Rm * S'(P X) * y
%    = Rm * S(y) * (X' kron eye(3)) p
%
% P.P = estimated projection matrix with uncertainty, 
%       using algebraic minimization, thus neglecting the accuracy
% P.Crr = CovM of reduced parameters, derived by variance propagation
%
% points may not be conditioned.
%
% Wolfgang Förstner 2/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_estimation_ml_Projection_3D_2D_from_point_pairs

function P = sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs(X,y)

% number of point pairs
N = size(X.e,1);
U = 12; 

% condition
[Xs,MX] = sugr_condition_Points(X);
[ye,My] = sugr_condition_Points(y);


% homogeneous coordinates
X = Xs.h;
% nonhomogeneous coordinates
y = ye.e;


%% estimate P algebraically
%
% build coefficient matrix
A = zeros(2*N,U);
 for n = 1:N
%     % use first two coordinates
    % S * Kronecker 3xU -matrix
    SK  =  calc_S([y(n,:),1]') * ...
            [X(n,1)*eye(3), X(n,2)*eye(3) X(n,3)*eye(3) X(n,4)*eye(3)];
    % 2xU part of A
    A(2*n-1:2*n,:)  = SK(1:2,:); 
 end

% partition
[Um,Dm,Vm] = svd(A,'econ');
log(abs(diag(Dm)));
% find algebraically best solution
h          = Vm(:,U);
% determine A+ = pseudo inverse of A with last singular value = 0
Di         = inv(Dm+eps*eye(U));
Di(U,U)    = 0;
A_plus     = Vm * Di * Um';

% rehape and choose unique sign
Pe = reshape(h,3,4);
Pe = Pe * sign(det(Pe(1:3,1:3)));

%% 
% determine covariance matrix Chh = A^+ Cgg A^+'

% determine Cgg of residuals of constraints from 
% Rm S(yh) P X = g                    yh = [y;1]
% d(Rm S(yh) P X) = d(g)
% Rm S(yh) P dX - Rm S(P X) dyh = dg
% Rm S(yh) P Jr dXr - Rm S(P X) Rm' dy = dg
%
CXyXy = sparse(zeros(2*N,2*N));
for n=1:N
    % VarProp for g 
    % [Ssy,Ry] = calc_S_reduced([y(n,:),1]');
%     S   = calc_S([y(n,:),1]');
    Rm  = [eye(2) zeros(2,1)];
    J_X = + Rm * calc_S([y(n,:),1]') * Pe * null(X(n,:));
    J_y = - Rm * calc_S(Pe * X(n,:)') *  Rm' ;
    % effect of both
    CXyXy(2*n-1:2*n,2*n-1:2*n)  = CXyXy(2*n-1:2*n,2*n-1:2*n) ...
        + J_y * squeeze(ye.Cee(n,:,:)) * J_y' ...
        + J_X * squeeze(Xs.Crr(n,:,:)) * J_X';                               %#ok<SPRIX>
end
Chh = A_plus * CXyXy * A_plus';

Pc = sugr_Projection_3D_2D(Pe,Chh);

% uncondition projection matrix
P = sugr_uncondition_Projection(Pc,My,MX);
return
