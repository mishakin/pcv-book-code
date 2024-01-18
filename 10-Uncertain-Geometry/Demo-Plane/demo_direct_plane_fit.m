%% demo plane fitting direct method
%
% generate points in [-1,+1]x[-0.5,+0.5] in arbitry plane
% add noise
%
% estimate plane parameters of centroid representation algebraically
%
% check estimation
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 06/18
% wenzel@igg.uni-bonn.de

% clear all 
close all
clearvars

addpath(genpath('../../General-Functions'))

disp('--------------------------------------------')
disp('Check direct estimation of plane from points')
disp('--------------------------------------------')

% number of samples N_samples
N_samples = 10000;
disp(['Number of samples                      : ', num2str(N_samples)])

% standard deviation of points
sigma_x = 1e-3;
disp(['MInimal standard deviation of points   : ', num2str(sigma_x)])

% density
pointSpacing = 0.1;
disp(['Point spacing in [-1,+1]x[-0.5,+0.5]   : ', num2str(pointSpacing)])
numPoints = length(-1:pointSpacing:1)*length(-0.5:pointSpacing:0.5);
disp(['Number of points                       : ', num2str(numPoints)])

%% random seed
%init_rand      = 1; 
init_rand      = 0;
init_rand_seed(init_rand);

%% additional control variables

% Parameterize an initial plane (z=1 plane) == true plane
iniPlane = [0, 0, 1, -1]'; % z=1 plane

% uncertainty
% random variances between 0.0001 and 0.001 

SIGMA = sigma_x * randi(10, [1, numPoints]);
%SIGMA = sigma_x * randInt(10, [1, numPoints]);
w     = 1./SIGMA.^2;

% alpha for checks
alpha = 0.001;

%% generate regularly spaced points on the initial plane == true points
[m1, m2] = meshgrid(-1:pointSpacing:1, -0.5:pointSpacing:0.5);
numPoints = numel(m1(:));
iniPoints = [m1(:), m2(:), ones(numPoints,1)]'; % points of z=1 plane

%% make up a random transformation == true transformation
randRot = pi/4 * (rand(1,3) - 0.5); % random rotation within [-90 90] degrees
randR = calc_Rot_r([randRot(1), randRot(2), randRot(3)]);
randT = 20 * (rand(3,1) - 0.5); % random translation within range [-10 10]

%% transform initial Plane == true plane
plane_nd = [randR, randT; 0 0 0 1]' * iniPlane; 
plane_nd = plane_nd * sign(-plane_nd(4)); % d is always positive

%% Sample data and store results
Omega_s     = zeros(1,N_samples);
d_prime_s   = zeros(1,N_samples);
normal_XY_s = zeros(2,N_samples);
var_q_s     = zeros(1,N_samples);
var_phi_s   = zeros(1,N_samples);
var_psi_s   = zeros(1,N_samples);
d_primes_s  = zeros(1,N_samples);

start = cputime;
for n_samples = 1:N_samples
    %transform points on the initial plane and add noise
    points = randR'*(iniPoints - repmat(randT, 1, numPoints)) + ...
        repmat(SIGMA, 3, 1) .* randn(3, numPoints);
    
    % fit a plane using direct method
    [Xo, Q, var_q, var_phi, var_psi, est_sigma_0_2] = ...
        sugr_direct_fit_plane_centroid_form_to_points(points, SIGMA.^2);
    
    Omega_s(n_samples)       = est_sigma_0_2 * (numPoints-3);
    d_primes_s(n_samples)    = -plane_nd(4) - plane_nd(1:3)' * Xo;
    normal_XY_s(:,n_samples) = Q(:,1:2)' * plane_nd(1:3);
    var_q_s(n_samples)       = var_q;
    var_phi_s(n_samples)     = var_phi;
    var_psi_s(n_samples)     = var_psi;
    
end
disp(['CPU time for ', num2str(N_samples),' samples             : ',num2str(cputime-start)]);
%% check
Red  = numPoints-3;

%% Check estimation

% collect parameters
est_par      = [d_primes_s;normal_XY_s];
% collect variances
CovM_true    = diag([var_q,var_phi,var_psi]);
%
check_estimation_result(Red,zeros(3,1),CovM_true,Omega_s/Red,est_par',1-alpha,' algebraic plane fit');






