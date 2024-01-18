%% Algebraic estimate of intersection Point_2D from N Lines_2D
%
% [x_est,sigma_0_est_D,sigma_0_est_O] =
%         sugr_algebraic_estimation_Point_2D_from_lines(lines,spherically)
%
% lines = struct
%         lines.h    = N x 3 matrix of sherically normalized homogeneous coordinates
%         lines.Crr  = N x 2 x 2 array of reduced covariance matrices
%         lines.type = N x 1 vector of 2's
% spherically = 1 use spherically normalized lines
%             = 0 use Euclideanly normalized lines
%
% x_est       = structure of algebarically estimated intersection point |x_est|=1
% sigma_0_est_O = estimated variance factor sqrt(v' Sigma^{-1} v /(N-2))
% sigma_0_est_D = estimated variance factor sqrt(D33^2 /(N-2))
%                 from SVD of points
%
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_Line_2D,  sugr_estimation_ml_Point_2D_from_Lines
% sugr_estimation_geometric_Point_2D_from_Lines


function [ex,sigma_0_est_D,sigma_0_est_O] = sugr_estimation_algebraic_Point_2D_from_Lines(lines,sph)

lh   = lines.h;
lCrr = lines.Crr;

[N,~] = size(lh);
Red = N-2;

% algebraic estimation
% Euclidean normalization changes
if sph==0
    lh=lh./(sqrt(lh(:,1).^2+lh(:,2).^2)*[1,1,1]);
end
[U,D,V]=svd(lh,'econ');
est_x = V(:,3);

% Covariance matrix: from c+v = A x, A = (l_n'), c = (c_n) = (l_n' x)
% dx = (A' A)^-1 A' c = A^+ c --> Cxx = A^+ Ccc A^+', Ccc = (x' Clnn x)

% Pseudoinverse of D
invD = zeros(3);
for i=1:2
    invD(i,i) =1/D(i,i);
end
%pseudoinvers of A
Ap = V * invD * U';                              % pseudo inverse (N x 3)

Omega = 0;
Cxx   = zeros(3);
for n=1:N
    Jln      = null(lh(n,:));                       % Jacobian l -> l_r
    Clnn     = Jln * squeeze(lCrr(n,:,:)) * Jln';    % Cl_nl_h
    cgn      = est_x' * lh(n,:)';                    % cgn
    var_cgn  = est_x' * Clnn * est_x;               % var(cgn)
%     c_sc =[cgn/sqrt(var_cgn)];
    Omega    = Omega + cgn^2/var_cgn;               %
    Cxx      = Cxx + var_cgn * Ap(:,n) * Ap(:,n)';  % Cxx update
end
if Red > 0
    sigma_0_est_D = sqrt(D(3,3)^2/Red);
    sigma_0_est_O = sqrt(Omega/Red);
else
    sigma_0_est_D = 1;
    sigma_0_est_O = 1;
end

ex = sugr_Point_2D(est_x,Cxx);

