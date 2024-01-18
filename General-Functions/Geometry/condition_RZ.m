%% condition Pose RZ
%
% Usage:
%   RZ = unconditioned Pose
%        RZ.RZ = [R,Z]
%        RZ.C  = D([dr;dZ])
%  
%   M    - 4x4 condition matrix scene points
% 
%   RZc  - Pose, assuming Pc=[I3,0]*M(Rc,Zc)
%        RZc.RZ = [Rc,Zc]
%        RZc.C  = D([drc;dZc])
%
%
% Wolfgang Förstner 1/2018
% wfoerstn@uni-bonn.de 
% 
% conditioning-matrix: Xec = (Xe - t)/s -> Xc = M * X 
% [Xch] = [eye(3)/s -t/s]   [Xh] = [Xh/s-t*X0/s]
% [Xc0]   [0'         1 ]   [X0]   [   X0      ]
% y = N(R*(X0/Xh-Z)) = N(R*((X0/Xh-t)/s -(Z-t)/s) = N(R*(Xc-Zc))
% Zc = (Z-t)/s 


function RZc = condition_RZ(RZ, M)

s = 1/M(1,1);
t = -M(1:3,4)*s;
% Conditioning of Z
RZc.RZ = [RZ.RZ(:,1:3), (RZ.RZ(:,4)-t)/s];

% variance propagation
Ji = [eye(3) zeros(3); zeros(3) eye(3)/s];
RZc.C = Ji*RZ.C*Ji';

return


