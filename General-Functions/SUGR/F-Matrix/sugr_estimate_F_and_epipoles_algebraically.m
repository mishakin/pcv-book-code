%% Algebraic estimate F-matrix and epipoles: from n > 7 points
%
% [F,Cff,el,Cll,er,Crr]=sugr_estimate_F_and_epipoles_algebraically(X,Cxx)
%
% X   = n x 4 Matrix of image points, not necessarily conditioned
% Cxx = n x 16 matrix of vectorized covariances of Euclidean point
%       coordinates, correspoinding points may be correlated
%       Cxx(i,:)= vec ( C(x_i'x_i')  C(x_i'x_i'') )
%                      C(x_i''x_i')  C(x_i''x_i'') )
%
% F   = 3x3 matrix
% Cff = 9x9 Cov(vec(F))
% el  = left epipole
% Cll = 3x3 CovM of el
% er  = right epipole
% Crr = 3x3 CovM of er
%
% Wolfgang Förstner 2/2012
% wfoerstn@uni-bonn.de

function [F, Cff, el, Cll, er, Crr] = sugr_estimate_F_and_epipoles_algebraically(X, Cxx)

% Number of observations N
[N, ~] = size(X);

% condition
% left
mxl = mean(X(:, 1)); sxl = std(X(:, 1));
myl = mean(X(:, 2)); syl = std(X(:, 2));
Y(:, 1) = (X(:, 1) - mxl) / sxl;
Y(:, 2) = (X(:, 2) - myl) / syl;
Tl = [1/sxl, 0, -mxl/sxl; 0, 1/syl, -myl/syl; 0 0 1];
% right
mxr = mean(X(:, 3)); sxr = std(X(:, 3));
myr = mean(X(:, 4)); syr = std(X(:, 4));
Y(:, 3) = (X(:, 3) - mxr) / sxr;
Y(:, 4) = (X(:, 4) - myr) / syr;
Tr = [1/sxr, 0, -mxr/sxr; 0, 1/syr, -myr/syr; 0, 0, 1];


% Coefficients for f=vecF
A = [Y(:, 1) .* Y(:, 3), ...
     Y(:, 2) .* Y(:, 3), ...
     1 * Y(:, 3), ...
     Y(:, 1) .* Y(:, 4), ...
     Y(:, 2) .* Y(:, 4), ...
     Y(:, 4), ...
     Y(:, 1), ...
     Y(:, 2), ...
     ones(N, 1)...
     ];
% best f, approximate
[U, D, V] = svd(A);

% Pseudoinverse for Covariance matrix
Dinv = D';
for i = 1:8
    Dinv(i, i) = 1 / D(i, i);
end
Ainv = V * Dinv * U';

% approximate F
Fa = reshape(V(:, 9), 3, 3);

% partition
[U, D, V] = svd(Fa);

% enforce singularity
D(3, 3) = 0;

% final estimate conditioned
F = U * D * V';
normF = norm(F, 'fro');
F = F/normF; % normalize such that Frobenius norm is 1
f = F(:);

% check
% diag([Y(:, 1), Y(:, 2), ones(N, 1)] * F * [Y(:, 3), Y(:, 4), ones(N, 1)]')

% Covariance matrix: Ainv * Cn * Ainv' restricted to det F=0
Cff = zeros(9);
var_n = zeros(N,1);
for n = 1:N
    % n-th constraint vector
    xx = [[Y(n, 3:4), 1] * F',[Y(n,1:2),1]*F];
    % n-th covariance matrix Euclidean
    C = reshape(Cxx(n, :), 4, 4);
    % n-th covarinace matrix homogeneous
    Cn = [C(1:2, 1:2) zeros(2, 1) C(1:2, 3:4) zeros(2, 1); ...
          zeros(1, 6); ...
          C(3:4, 1:2) zeros(2, 1) C(3:4, 3:4) zeros(2, 1); ...
          zeros(1, 6)...
          ];
    % Jacobian for x->y
    J = [Tl, zeros(3); zeros(3), Tr];
    % covariance matrix for conditioned homogeneous coordinates
    Cn = J * Cn * J';
    % variance of constraint
    var_n(n) = (xx * Cn * xx');
    % covarance matrix of vec f approximate
    Cff = Cff + Ainv(:, n) * Ainv(:, n)' * var_n(n);
end
% rankCffa = rank(Cff)

% projector for enforcing det and norm constraint
FA = adjunctMatrix(F)';
H = [FA(:), f];
P = eye(9) - H * inv(H'*H)*H';                                             %#ok<MINV>
% covariance matrix for conditioned F
Cff = P * Cff * P;

% rankCff = rank(Cff)

% uncoditioned F
F = Tl' * F * Tr;

% covariance matrix of unconditioned F
J = kron(Tr',Tl');
Cff = J * Cff * J';

% rank(Cff);

% determine epipoles and their covariances
% left
el = cross(F(:, 1), F(:, 2));
fac = norm(el);
el = el / fac;
Jl = [-calc_S(F(:, 2)), calc_S(F(:, 1))];
Cll = Jl * Cff(1:6, 1:6) * Jl'/fac^2;
%right
er = cross(F(1, :), F(2, :))';
fac = norm(er);
er = er / fac;
Jr = [-calc_S(F(2, :)), calc_S(F(1, :))];
Crr = Jr * Cff([1, 4, 7, 2, 5, 8], [1, 4, 7, 2, 5, 8]) * Jr'/fac^2;

% check
%[X(:,1),X(:,2),ones(N,1)]*F*[X(:,3),X(:,4),ones(N,1)]'

