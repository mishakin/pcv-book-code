%% Algebraic estimation of 2D homography from point pairs, 
% assuming independent observations
% and conditioned points
%
% H = sugr_estimation_algebraic_Homography_2D_from_point_pairs(P) 
%
% P contains point pairs
% P.h = [X,Y]
% * X = N x 3 matrix if points 
% * Y = N x 3 matrix if points
% P.Crr = N x 2 x 2 contains reduced covariance matrix
% 
%
% 0 != S^s(y) * H * x 
%    = Rm(y) * S(y) * H * x 
%    = Rm(y) * S'(H x) * y
%    = Rm(y) * S(y) * (x' kron eye(3)) h 
%
% H.H   = estimated homography with uncertainty, 
%         using algebraic minimization, thus neglecting the accuracy
% H.Crr = CovM of reduced parameters, derived by variance propagation
%
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% wf 2/2010

function H = sugr_estimation_algebraic_Homography_2D_from_point_pairs(P)

% pair coordinates
Ph   = P.h;
% pair CovM reduced
PCrr = P.Crr;
% homogeneous coordinates
X = Ph(:,1:3);
Y = Ph(:,4:6);
% number of point pairs
N = size(X,1);

%% 
% estimate H algebraically
%
% build coefficient matrix
A = zeros(2*N,9);
for n = 1:N
    % select 2x3 skew matrix
    [Ssy,~] = calc_S_reduced(Y(n,:)');
    % Kronecker 3x9 -matrix
    K  = [X(n,1)*eye(3), X(n,2)*eye(3) X(n,3)*eye(3)];
    % 2x9 part of A
    A(2*n-1:2*n,:)  = Ssy * K; 
end

% partition
N_rows = size(A,1);
if N_rows == 8
    h = null(A);
    A_plus = A'/(A*A');
else
    [U,D,V] = svd(A,'econ');
    % find algebraically best solution
    h        = V(:,9);
    % determine A+ = pseudo inverse of A with last singular value = 0
    Di       = inv(D+eps*eye(9));
    Di(9,9)  = 0;
    A_plus   = V * Di * U';
end
% choose unique sign (largest element > 0)
[~,imax] = max(abs(h));
h        = h * sign(h(imax));
He       = reshape(h,3,3);
He=He/(abs(det(He))^(1/3));


%% determine covariance matrix Chh = A^+ Cgg A^+'

% determine Cgg of residuals g of constraints from 
% S^s(y) H x                        = g 
% Ry Sy  H x                        = g
% d(Ry Sy H x)                      = d(g)
% Ry Sy H dx     - Ry S(H x) dy     = dg
% S^sy  H Jr dxr - Ry S(H x) Jr dyr = dg
%

%Cxyxy = sparse(zeros(2*N,2*N));
Cxyxy = spalloc(2*N,2*N,2*N);
for n=1:N
    % VarProp for g 
    [Ssy,Ry] = calc_S_reduced(Y(n,:));
    J_x = + Ssy * He * null(X(n,:));
    J_y = - Ry * calc_S(He * X(n,:)') *  null(Y(n,:));
    % effect of both
    J = [J_x J_y];  % 2 x 4 Jacobian res/pr = [res/xr res/yr]
    Cxyxy(2*n-1:2*n,2*n-1:2*n)  = Cxyxy(2*n-1:2*n,2*n-1:2*n) ...
        + J * squeeze(PCrr(n,:,:)) * J';                                    %#ok<SPRIX>
end
Chh = A_plus * Cxyxy * A_plus';

H = sugr_Homography_2D(He,Chh);
