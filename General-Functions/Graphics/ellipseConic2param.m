%% estimates values of parametric represenation of a 2D ellipse
%
% Usage:
%   p = ellipseConic2param(C)
%
%   C: double 3x3 ellipse-conic, det(Chh)!>0
%   p: double 5x1, [x0,y0,a,b,phi]
%
% 10.12.13 Susanne Wenzel
% wenzel@igg.uni-bonn.de

function p = ellipseConic2param(C)

% help values
Chh = C(1:2,1:2);

if det(Chh)<=0
    keyboard
end

Chhinv = inv(Chh);
ch0 = C(1:2,  3);
c00 = C(  3,  3);

% centre point
x0 = -Chhinv*ch0;

% eigenvalue decomposition of normalized homogeneous part of the conic
Chhe = -Chh/(c00-ch0'*Chhinv*ch0); %#ok<MINV>
[R,Lambda] = eig( Chhe );

l = diag(Lambda);
[minlambda,idx_min] = min(l);
maxlambda = max(l);

% main radii
a = sqrt(1/minlambda);
b = sqrt(1/maxlambda);

if  abs(minlambda-maxlambda)<eps
    % its a circle, orientation 0
    phi = 0;
else
    phi = atan2(R(2,idx_min),R(1,idx_min));
    if phi<0
        phi = phi+2*pi;
    end
end

 p = [x0',a,b,phi];