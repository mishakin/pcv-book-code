%% determines constraint value and Jacobians for 2D line estimation
%
% g= ^l2 - a ^l1 -b = 0  (== y_i-ax_i-b=0,i=1,...,N)
%
% lv     = 2*G vector of observations [x;y]
% lnu    = 2*G vector of approximate fitted observations [x^a;y^a]
% xa     = 2x1 vector, approximate values for parameters xa=[a;b]
%
% cg     = Gx1 vector of residual of constraints
% A      = GxU Jacobian dg/dx
% B      = NxG Jacobian (dg/dl)'
%
% Wolfgang Förstner
% wfoerstn@uni-bonn.de

function  [cg,A,B] = cg_2D_line(lv,lnu,xa)

% number of constraints
G2 = size(lv,1);
G  = G2/2;                        % == N

% Jacobians
A  = - [lv(1:G) ones(G,1)];       % G x U
B = [-xa(1)*eye(G) eye(G)]';      % N x G

% residual of constraints
cg = -(lv(G+1:2*G) - xa(1) * lv(1:G) - xa(2)) + B'*(lnu-lv);  % G x 1


