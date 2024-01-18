%% Fits a plane in centroid form to a set of 3D points
%
% [Xo, Q, var_q, var_phi, var_psi] = sugr_direct_fit_plane_centroid_form_to_points(X, CovX)
%
% * X    = 3 x n matrix, n 3D points that lie on a plane
% * CovX = 1 x n or 3 x n or 9 x n matrix, covariance matrices of the points.
% Each column of CovX is either the variance of a point, or the diagonal/all elements of the covariance matrix of a point.
%
% * Xo  = 3 x 1 vector, centroid of the fitted plane
% * Q   = 3 x 3 matrix, rotation matrix of the fitted plane. N = Q(:,:,3);
% * var_q   = 1 x 1 scalar, theoretical variance of fitted plane along the normal
% * var_phi = 1 x 1 scalar, theoretical variance of X component of the normal
% * var_psi = 1 x 1 scalar, theoretical variance of Y component of the normal
% * est_sigma_0_2   = 1 x 1 scalar, estimated variance factor
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% change log
% 29/08/2016    wf: added est_sigma_0_2 as output
% 26/08/2016    first version
% 28/08/2016    CovX can now be a 1 x n, 3 x n or 9 x n matrix

function [Xo, Q, var_q, var_phi, var_psi, est_sigma_0_2] = sugr_direct_fit_plane_centroid_form_to_points(X, CovX)

%% check input

% check that we dont have too few arguments
if nargin < 1, error('Not enough input arguments. '), end

% or too many ...
if nargin > 2, error('Too many input arguments'), end

% otherwise we have at least points; check if they are 3D
[dim, numPoints] = size(X);
if dim ~= 3, error('Input should be a 3 x n matrix containing 3D points. '), end

% check if we have enough points
if numPoints < 3, error('At least three points are needed for plane fitting. '), end

% if variances are not given use default value of 1 for all points
if nargin == 1
    warning('A default variance of 1 will be used for all points. ')
    CovX = repmat(ones(3,1), 1, numPoints);
end

% check the form of (co)variances if given
if nargin == 2
    
    [dim2, np2] = size(CovX);
    
    % check that variances are given for all points
    if np2 ~= numPoints, error('Number of variances does not match the number of points. '), end
    
    % if for each point only one variance value is given then all is fine
    if dim2 == 1
        CovX = repmat(CovX, 3, 1);
        
        % if for each point three variances are given then issue a warning but we can proceed
    elseif dim2 == 3
        warning('For anisotropic covariance matrices the result may not be optimal. '),
        
        % if for each point all 9 elements of a covariance matrix are given issue a warning and discard the off-diagonal elements
    elseif dim2 == 9
        warning('Off-diagonal elements of covariance matrices will be ignored. '),
        %CovX = CovX([1 5 9], :);
        CovX = repmat((CovX(1, :)+CovX(5, :)+CovX(9, :))/3,3,1);
        
        % hopefully one should never end up here!
    else
        error('Input variances/covariances should be in the form of a 1 x n, 3 x n or 9 x n matrix. '),
    end
end

%% the weights
w = 1 ./ CovX;

%% the weighted centroid
Xo = sum(w .* X, 2) ./ sum(w, 2);

%% the moment matrix
M = w .* (X - repmat(Xo, 1, numPoints)) * (X - repmat(Xo, 1, numPoints))' ;

%% the rotation matrix Q
[Q, G] = eig(M);

% sort eigenvalues in descending order
[srtval, srtid] = sort(diag(G), 'descend');
G = diag(srtval);

% sort eigenvectors correspondingly
Q = Q(:, srtid);

% if Q is a reflection convert it to a rotation
Q = Q*diag([1,1,det(Q)]);

%% the variances of plane parameters
w_bar = mean(w(:));

% redundancy
Red = numPoints - 3;

% estimated Omega
Omega = G(3,3);

if Red > 0
    est_sigma_0_2 = Omega/Red;
else
    est_sigma_0_2 = 1;
end

% theoretical variance along normal
var_q   = 1 / numPoints / w_bar ;

% theoretcial variance of x component of the normal
var_phi = 1 / G(1, 1)  ;

% theoretical variance of y component of the normal
var_psi = 1 / G(2, 2)  ;


