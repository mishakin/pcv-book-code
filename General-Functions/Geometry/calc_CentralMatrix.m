%% CentralMatrix: 
% determines Jacobian dm/dx for determining mean 
% of non-homogeneous vectors x and y 
% 
% Usage:
% Z = calc_CentralMatrix(x)
%
%   x - D-vector
% 
%   Z - CentralMatrix(x), such that 
%       Central matrix of d-vector Midpoint z = Z(x) y of x, y
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de 

function Z = calc_CentralMatrix(x)
    D = length(x);
    Z = x(D) * eye(D) + x * [zeros(1,D-1) 1];
