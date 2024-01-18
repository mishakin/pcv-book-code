%% generate_observed_ARp: generates samples of autoregressive processes
%
% [x,y,sigma_x,sigma_y] = 
%   generate_observed_ARp(N,p,q,sigma_e,sigma_n);
% 
% N         = number of points
% p         = order of process
% q         = decay (abs(q) < 1)
% sigma_e   = process noise
% sigma_n   = observation noise
%
% x         = true sample
% y         = noisy sample
% sigma_x   = theoretical std. of x
% sigma_y   = theoretical std. of y
%
% Wolfgang Förstner 01/2015
% wfoerstn@uni-bonn.de

function [x,y,sigma_x,sigma_y] = generate_observed_ARp(N,p,q,sigma_e,sigma_n)

% initiate sample
x = zeros(N,1);
y = zeros(N,1);

% determine coefficients from polynomial
c = poly(q*ones(1,p));        % p roots q -> polynomial coeff.
a = -c(2:p+1)';               % coeff. of ARp = polynomial coeff.

% theoretical std of x and y
sigma_x = sigma_e/sqrt((1-q)^p);
sigma_y = sqrt(sigma_x^2+sigma_n^2);

% generate true and noisy sample
for n=p+1:N
    x(n) = sum(a(1:p).*x(n-1:-1:n-p))+randn(1)*sigma_e;
    y(n) = x(n)+randn(1)*sigma_n;
end


