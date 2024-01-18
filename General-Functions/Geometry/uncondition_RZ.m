%% uncondition Pose RZ
%
% Usage:
%   RZc = unconditioned Pose
%        RZc.RZ = [Rc,Zc]
%        RZc.C  = D([drc;dZc])
%  
%   M   - 4x4 condition matrix from scene points
% 
%   RZ  - Pose, assuming P=[I3,0]*M(R,Z)
%        RZ.RZ = [R,Z]
%        RZ.C  = D([dr;dZ])
%
% Wolfgang Förstner 1/2018
% wfoerstn@uni-bonn.de 
%
% conditioning-matrix: Xec = (Xe - t)/s -> Xc = M * X 
% [Xch] = [eye(3)/s -t/s]   [Xh] = [Xh/s-t*X0/s]
% [Xc0]   [0'         1 ]   [X0]   [   X0      ]
% y = N(R*(X0/Xh-Z)) = N(R*((X0/Xh-t) -(Z-t))/s) = N(R*(Xc-Zc))
% Zc = (Z-t)/s  -> Z = s*Zc+t
% 

function RZ = uncondition_RZ(RZc, M)
% Parameters of conditioning
s =  1/M(1,1);
t = -M(1:3,4)*s;

% Unonditioning of Z
RZ.RZ = [RZc.RZ(:,1:3), s*RZc.RZ(:,4)+t];

% variance propagation
J = [eye(3) zeros(3); zeros(3) eye(3)*s];
RZ.C = J*RZc.C*J';


return

