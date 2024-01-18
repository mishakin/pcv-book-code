%% Test singularity of line through four lines
% Given four lines L_j,j=1,2,3,4, determine lines M_i passing these
%
% See PCV Section 7.2.1.3 Three and More Entities
%
% 0 or 2 solutions
%
% Cases:
%   random
%   three parallel, one other line
%   three intersecting lines, one other
%   three coplanar intersecting lines, one other
%   two intersecting lines, two others
%
% Output per case
%   3D lines
%   nullspace
%   coefficients (should not be zero) see PCV (7.53)
%   solutions (may be identical)
%   two lines, if possible
%
% Wolfgang Förstner 10/2016
% wfoerstn@uni-bonn.de
%
% last changes: Susanne Wenzel 12/17
% wenzel@igg.uni-bonn.de

% clear all
clearvars
clc
close all

addpath(genpath('../General-Functions'))

fprintf('\n------ DEMO Line Through Four Lines ------\n')

fprintf('\nIn order to shorten the output, we display transposed vectors\n')


% four random lines
disp('------------- four random lines --------------')
L1 = calc_Pi(rand(4,1))*randn(4,1);
L2 = calc_Pi(rand(4,1))*randn(4,1);
L3 = calc_Pi(rand(4,1))*randn(4,1);
L4 = calc_Pi(rand(4,1))*randn(4,1);
M = [L1,L2,L3,L4]                                                           %#ok<NOPTS>
disp(['rank(M) = ', num2str(rank(M))])

% determine polynomial for intersecting line
fprintf('\ndetermine polynomial for intersecting line\n')
NM = null(M'*calc_Dual);
null_space_M = NM                                                           %#ok<NOPTS>

N1 = NM(:,1);
N2 = NM(:,2);
a = (N1-N2)'*calc_Dual*(N1-N2);
b = 2*(N1-N2)'*calc_Dual*N2;
c = N2'*calc_Dual*N2;
coefficients = [a,b,c]                                                       %#ok<NOPTS>
rr = roots([a,b,c]);
disp(['solutions: ',num2str(rr')])

M1 = rr(1)*N1+(1-rr(1))*N2;
M2 = rr(2)*N1+(1-rr(2))*N2;
two_lines = [M1';M2'];
fprintf('\ntwo lines \n')
disp(['M_1 = [',num2str(M1'),']'])
disp(['M_2 = [',num2str(M2'),']'])
incidences = two_lines*calc_Dual*M                                          %#ok<NASGU,NOPTS>

% three parallel
fprintf('\n------------- three random but parallel lines, one other --------------\n')
% randomly generate two directions
D123 = [rand(3,1);0];
D4 = [rand(3,1);0];
% generate 3 3D lines having the same normal direction
L1 = calc_Pi(rand(4,1))*D123;
L2 = calc_Pi(rand(4,1))*D123;
L3 = calc_Pi(rand(4,1))*D123;
% generate another 3D line with different normal direction
L4 = calc_Pi(rand(4,1))*D4;
L1 = L1/L1(1);
L2 = L2/L2(1);
L3 = L3/L3(1);
L4 = L4/L4(1);
M = [L1,L2,L3,L4]                                                           %#ok<NOPTS>
disp(['rank(M) = ', num2str(rank(M))])

NM = null(M'*calc_Dual);
null_space = NM                                                             %#ok<NASGU,NOPTS>
N1 = NM(:,1);
N2 = NM(:,2);
a = (N1-N2)'*calc_Dual*(N1-N2);
b = 2*(N1-N2)'*calc_Dual*N2;
c = N2'*calc_Dual*N2;
coefficients_polynomial_no_solution = [a,b,c]                               %#ok<NASGU,NOPTS>
solutions = roots([a,b,c])';
disp(['solutions: ',num2str(solutions)])

% three through one point
fprintf('\n------------- three random but intersecting lines, one other --------------\n')
% define 2 3D points
X1 = rand(4,1);
X2 = rand(4,1);

% generate random 3d lines passing through the same point X1
L1 = calc_Pi(rand(4,1))*X1;
L2 = calc_Pi(rand(4,1))*X1;
L3 = calc_Pi(rand(4,1))*X1;
% generate another random 3d lines passing through X2
L4 = calc_Pi(rand(4,1))*X2;
M = [L1,L2,L3,L4]                                                           %#ok<NOPTS>
disp(['rank(M) = ', num2str(rank(M))])

NM = null(M'*calc_Dual);
nullspace = NM                                                              %#ok<NOPTS,NASGU>
N1 = NM(:,1);
N2 = NM(:,2);
a = (N1-N2)'*calc_Dual*(N1-N2);
b = 2*(N1-N2)'*calc_Dual*N2;
c = N2'*calc_Dual*N2;
coefficients_polynomial_no_solution = [a,b,c]                               %#ok<NOPTS,NASGU>
solutions = roots([a,b,c])';
disp(['solutions: ',num2str(solutions)])

%% three coplanar three through one point
fprintf('\n------------- three coplanar intersecting lines, one other --------------\n')
% define three points on the plane
X = [1,0,0,1]';
Y = [0,1,0,1]';
Z = [0,0,1,1]';
% define another arbitrary point
T = [randn(1,3),1]';
% check whether they occasionally lie on a plane
disp(['Determinant should not be zero: ', num2str(det([X,Y,Z,T]))])

% 3D lines trough points on plane (X,Y,Z)
L1 = calc_Pi(X)*Y;
L2 = calc_Pi(Y)*Z;
L3 = calc_Pi(Z)*X;
% 3D line not on the same plane
L4 = calc_Pi(X)*T;

M = [L1,L2,L3,L4]                                                           %#ok<NOPTS>
disp(['rank(M) = ', num2str(rank(M))])

NM = null(M'*calc_Dual);
nullspace = NM                                                              %#ok<NOPTS>
N1 = NM(:,1);
N2 = NM(:,2);
as = (N1-N2)'*calc_Dual*(N1-N2);
bs = 2*(N1-N2)'*calc_Dual*N2;
cs = N2'*calc_Dual*N2;
coefficients_polynomial_no_solution = [as,bs,cs]                            %#ok<NOPTS>
solutions = roots([as,bs,cs])';
disp(['solutions: ',num2str(solutions)])


%% two intersecting lines
fprintf('\n------------- two intersecting lines, two others --------------\n')
% define 4 3D points as points on target lines
X = [-1,0,2,1]';
Y = [1,0,3,1]';
Z = [0,1,2,1]';
T = [1,-1,3,1]';
% Line in direction of X-axis
L1 = [1,0,0,0,0,0]';
% Line in direction of Y-axis
L2 = [0,1,0,0,0,0]';
% Line through points X and Y, intersection the X-axis
L3 = calc_Pi(X)*Y;
% Line through points Z and T, intersecting the Y-axis
L4 = calc_Pi(Z)*T;

M = [L1,L2,L3,L4]                                                           %#ok<NOPTS>
disp(['rank(M) = ', num2str(rank(M))])

NM = null(M'*calc_Dual);
null_space = NM                                                             %#ok<NOPTS>
N1 = NM(:,1);
N2 = NM(:,2);
a = (N1-N2)'*calc_Dual*(N1-N2);
b = 2*(N1-N2)'*calc_Dual*N2;
c = N2'*calc_Dual*N2;
coefficients_polynomial = [a,b,c]                                           %#ok<NOPTS>
solutions = roots([a,b,c])';
disp(['solutions: ',num2str(solutions)])

M1 = rr(1)*N1+(1-rr(1))*N2;
M2 = rr(2)*N1+(1-rr(2))*N2;
two_lines = [M1';M2'];
fprintf('\ntwo lines \n')
disp(['M_1 = [',num2str(M1'),']'])
disp(['M_2 = [',num2str(M2'),']'])
incidences = two_lines*calc_Dual*M                                          %#ok<NOPTS>