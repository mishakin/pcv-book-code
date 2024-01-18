%% principal_distance_from_F: determina single c from fundamental matrix
%
% c = principal_distance_from_F(F)
%
% F = fundamental matrix, rank 2-matrix
%     assumed for conditioned cooridates
% 
% c = principal distance
%
% Reference:
% P. Sturm, Z.L. Cheng, P.C.Y. Chen, A.N. Poo 
% Focal length calibration from two views:
% method and analysis of singular cases
% Computer Vision and Image Understanding 99 (2005) 58–95
%
% Wolfgang Förstner 3/2013
% wfoerstn@uni-bonn.de


function c = principal_distance_from_F(F)

% enforce long principal distance as prior
factor = 1;
K0 = diag([factor,factor,1]);
F  = K0*F*K0;

% svd
[U,D,V] = svd(F);
a = D(1,1);
b = D(2,2);

%% build coefficients

% 1. linear equation
a1 = a*U(3,1)*U(3,2)*(1-V(3,1)^2)+b*V(3,1)*V(3,2)*(1-U(3,2)^2);
a0 = U(3,2)*V(3,1)*(a*U(3,1)*V(3,1)+b*U(3,2)*V(3,2));

% 2. linear equation
b1 = a*V(3,1)*V(3,2)*(1-U(3,1)^2)+b*U(3,1)*U(3,2)*(1-V(3,2)^2);
b0 = U(3,1)*V(3,2)*(a*U(3,1)*V(3,1)+b*U(3,2)*V(3,2));

% quadratic equation
c2 = a^2*(1-U(3,1)^2)*(1-V(3,1)^2)-b^2*(1-U(3,2)^2)*(1-V(3,2)^2);
c1 = a^2*(U(3,1)^2+V(3,1)^2-2*U(3,1)^2*V(3,1)^2)-...
     b^2*(U(3,2)^2+V(3,2)^2-2*U(3,2)^2*V(3,2)^2);
c0 = a^2*U(3,1)^2*V(3,1)^2-b^2*U(3,2)^2*V(3,2)^2;

f_squares = roots([c2,c1,c0]);
% fs = sqrt(f_squares);

% flin1 = -a0/a1;
% flin2 = -b0/b1;

c = 0;
for i = [1,2]
    f2 = f_squares(i);
    if isreal(f2) && f2 > 0
        if abs(a1*f2+a0) < 10^(-10) && abs(b1*f2+b0) < 10^(-10)
           c = sqrt(f2);
        end
    end
end

c = factor*c;


