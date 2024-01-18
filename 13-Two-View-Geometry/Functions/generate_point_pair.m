%% generate_point_pair
%
% X        = 4x1 vector, 3D point
% b,Rot    = relative orientation
% sigma    = radial error of direction
%
% [u,v] = generate_point_pair(X,b,R,sigma)
%
% Wolfgang Förstner 04/2014
% wfoerstn@uni-bonn.de


function [u,v] = generate_point_pair(X,b,Rot,sigma)

% projection matrices
P1 =     [eye(3) zeros(3,1)];
P2 = Rot*[eye(3),-b];

% true directions
u_true = P1*X; u_true = u_true/norm(u_true);
v_true = P2*X; v_true = v_true/norm(v_true);

% perturbed directions, only two angles are affected, 
R_1 = calc_Rot_rod(sigma*randn(3,1));
u   = R_1*u_true;
u   = u/norm(u);
R_2 = calc_Rot_rod(sigma*randn(3,1));
v   = R_2*v_true;
v   = v/norm(v);
