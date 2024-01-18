%% variance propagation by numerical differentiation
%
% Usage: 
%   [my,Syy,J]=var_prop_classical(@f,mx,Sxx,p)%
%
% Input:
%    f   = nonlinear function y=f(x,p) as handle of m-file
%    mx  = mean of x
%    Sxx = covaraince of x
%    p   = parameters for function
%
% Output 
%    my  = mean of y
%    Syy = covariance of y
%    J = Jacobian matrix
%
% condition: Sxx must not have zero diagonals, otherwise numerically unstable
%
function [my,Syy,J] = var_prop_classical(f,mx,Sxx,p)

% determine 
% mean of output
% dimension of output (must be coded in function)
if nargin == 3
    my = f(mx);
else
    my = f(mx,p);
end
nf = size(my,1);

% determine dimension of input
n = size(mx,1);

% factor for deviations from mean for Jacobian
k = 0.1;
 
% Jacobian
J = zeros(nf,n);
for ii = 1:n
    % finite difference of input
    s = sqrt(Sxx(ii,ii));
    % vector of effect: in direction if inout
    eff = zeros(n,1);
    eff(ii) = k;
    % symmetric finite difference
    if nargin == 3
        t = (f(mx+s*eff)-f(mx-s*eff))/(2*k*s+eps);
    else
        t = (f(mx+s*eff,p)-f(mx-s*eff,p))/(2*k*s+eps);
    end
    % Jacobian element
    J(:,ii) = t(:);
end;

% Variance
Syy = J*Sxx*J';
