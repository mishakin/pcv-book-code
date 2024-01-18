%% Covariance matrix of parameters of essential matrix
%
% CovM = sugr_CovM_E_matrix(X, Cpp, B, R)
%
% model
%        xl' * E * xr = 0
%        E = S(B) R'
%        xl = [I |  0] X
%        xr = [R | -B] X
%
% X   = N x 4 matrix of homogeneous 3D points, N >= 5
% Cpp = N x 4 x 4 Covariance matrix of reduced image coordinate pairs
% R   = 3 x 3 rotation matrix
% B   = 3 x 1 direction of base line, 1-->2
%
% Checked in sugr_estimation_ml_E_Matrix_from_point_pairs.m
%
%  see: PCV 13.3.5.1
%
% Wolfgang Förstner 2016-09-01
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

function CovM = sugr_CovM_E_matrix(X, Cpp, B, R)

% number of points
NN = size(X, 1);

% essential matrix
E = calc_S(B) * R';

% projection matrix right
Pr = R * [eye(3), - B];

n = 0;
% fitted image directions
for nn = 1:NN
    if norm(X(nn, :)) > 10 ^ (- 10)
        n = n + 1;
        xl(n, :) = X(nn, 1:3) / norm(X(nn, 1:3));                          %#ok<*AGROW>
        xr(n, :) = X(nn, :) * Pr';
        xr(n, :) = xr(n, :) / norm(xr(n, :));
    end
end
N = n;

% ancillary vectors
ll = xr * E';
lr = xl * E;
x1 = xr * R;

at = zeros(1, 5);
bt = zeros(1, 4);

% dtermine variances of constraints
btCb = zeros(N, 1);
for n = 1:N
    bt(1:2) = ll(n, :) * null(xl(n, :));
    bt(3:4) = lr(n, :) * null(xr(n, :));
    C = reshape(Cpp(n, :, :), 4, 4);
    btCb(n) = bt * C * bt';
end

% determine information matrix
Nm = zeros(5);
for n = 1:N
    at(1:2) = (calc_S(x1(n, :)') * xl(n,:)')' * null(B');
    at(3:5) = (calc_S(lr(n, :)') * xr(n,:)')';
    Nm = Nm + at' * at / btCb(n);
end

CovM = inv(Nm);



