%
% RRS_3POINT Direct solution for the spatial resection from three points
%
% [R,Z] = rrs_3point(X, m)
%
% Input:
%  X        - 3x3-matrix, rows are scene point coordinates
%  m        - 3x3-matrix, rows are directons in the camera system
%
% Output:
%  R      - Mx3x3-Matrix, M solution for the rotation from the scene into
%                         the camera system
%  Z      - Mx3-Vektor, M solutions for the projection centre 
%           We have m_i ~ R_j (X_i - Z_j)
%
%
% Autor: Christian Beder
%
function [R,Z] = rrs_3point(X, m)
% sides of scene triangle and angles between camera rays
a = norm(X(3,:) - X(2,:));
alpha = acos(m(3,:) * m(2,:)');
b = norm(X(1,:) - X(3,:));
beta = acos(m(1,:) * m(3,:)');
c = norm(X(2,:) - X(1,:));
gamma = acos(m(2,:) * m(1,:)');
% coefficients
A4 = ((a^2 - c^2) / b^2 - 1)^2 - 4*c^2/b^2 * cos(alpha)^2;
A3 = 4*((a^2-c^2)/b^2*(1-(a^2-c^2)/b^2)*cos(beta) - ...
        (1-(a^2+c^2)/b^2)*cos(alpha)*cos(gamma) + ...
        2*c^2/b^2*cos(alpha)^2*cos(beta));
A2 = 2*(((a^2-c^2)/b^2)^2 - 1 + 2*((a^2-c^2)/b^2)^2*cos(beta)^2 + ...
        2*(b^2-c^2)/b^2*cos(alpha)^2 - ...
        4*(a^2+c^2)/b^2*cos(alpha)*cos(beta)*cos(gamma) + ...
        2*(b^2-a^2)/b^2*cos(gamma)^2);
A1 = 4*(-(a^2-c^2)/b^2*(1+(a^2-c^2)/b^2)*cos(beta) + ...
        2*a^2/b^2*cos(gamma)^2*cos(beta) - ...
        (1-(a^2+c^2)/b^2)*cos(alpha)*cos(gamma));
A0 = (1+(a^2-c^2)/b^2)^2 - 4*a^2/b^2*cos(gamma)^2;
% roots of polynomial
v = roots ([A4, A3, A2, A1, A0]);
solutions=0;
% select real roots
for i=1:length(v)
    if isreal(v(i)) && (v(i) > 0)
        solutions=solutions+1;
        s(solutions,1) = sqrt(b^2 / (1+v(i)^2-2*v(i)*cos(beta)));
        s(solutions,3) = v(i)*s(solutions,1);
        xxx = s(solutions,3)*cos(alpha);
        yyy = sqrt(s(solutions,3)^2*cos(alpha)^2+a^2-s(solutions,3)^2);
        s(solutions,2) = xxx + yyy;
        if (xxx > yyy)
            solutions = solutions+1;
            s(solutions,[1,3]) = s(solutions-1,[1,3]);
            s(solutions,2) = xxx - yyy;
        end;
    end;
end;
% initiate absolute orientation
R = zeros(solutions, 3,3);
Z = zeros(solutions, 3);
b_o = X(3,:)' - X(1,:)';
c_o = X(2,:)' - X(1,:)';
R_o = [b_o / norm(b_o), cross(b_o,c_o)/norm(cross(b_o,c_o)), ...
        cross(b_o, cross(b_o, c_o)) / norm(cross(b_o, cross(b_o, c_o)))];
% for each solution determine R and Z
for i=1:solutions
    for j=1:3
        Xk(j,:) = s(i,j) * m(j,:);
    end;
    b_k = Xk(3,:)' - Xk(1,:)';
    c_k = Xk(2,:)' - Xk(1,:)';
    R_k = [b_k / norm(b_k), cross(b_k,c_k)/norm(cross(b_k,c_k)), ...
        cross(b_k, cross(b_k, c_k)) / norm(cross(b_k, cross(b_k, c_k)))];
    RR = R_k*R_o';
    R(i,:,:) = RR;
    Z(i,:) = (X(1,:)' - RR'*Xk(1,:)')';
end;
