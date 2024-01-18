%% Geometric estimate of intersection Point_2D from N Lines_2D
%
% [x_est, sigma_0, x] = sugr_estimation_geometric_Point_2D_from_lines(lines)
%
% lines = struct
%         lines.h    = N x 3 matrix of sherically normalized homogeneous coordinates
%         lines.Crr  = N x 2 x 2 array of reduced covariance matrices
%         lines.type = N x 1 vector of 2's
%
% x_est   = structure of geometricall estimated intersection point |x_est|=1 
% sigma_0 = estimated variance factor sqrt(v' Sigma^{-1} v /(N-2))
% x       = 2x1 estimated point
%
% Wolfgang Förstner 7/2012
% wfoerstn@uni-bonn.de
%
% See also sugr_Point_2D, sugr_Line_2D,  sugr_estimation_ml_Point_2D_from_Lines
% sugr_estimation_algebraic_Point_2D_from_Lines

function [xe,sigma_0,x] = sugr_estimation_geometric_Point_2D_from_Lines(lines)

% homogeneous coordinates
lh   = lines.h;
% number of edges
[N,~] = size(lh);
% redundancy
R = N-2;
% initiate normals
nv = zeros(2,1);
Nm = zeros(2);
di=zeros(N,2);
x0=zeros(N,2);
% convert to centroid from and build Normals
w = zeros(N,1);
for n=1:N
    sugr_l.h    = lh(n,:)';
    sugr_l.Crr  = reshape(lines.Crr(n,:,:),2,2);
    sugr_l.type = 2;
    [e,Cee] = sugr_Line_2D_hom2Hes(sugr_l);
    [x0n,pn,~,~]= sugr_Line_2D_Hes2cen(e,Cee);
   
    dirn =[cos(pn-pi/2),sin(pn-pi/2)]';
    wn   = 1;%1/sqn^2;
    w(n) = wn;
    Wn   = wn*(eye(2)-dirn*dirn');
    nv = nv + Wn*x0n;
    Nm = Nm + Wn;
    x0(n,:)=x0n';
    di(n,:)=dirn';
end
% covariance matrix
Cxx   = inv(Nm);
% estimated point
x     = Cxx * nv;                                                          %#ok<*MINV>
% sum of weighted squared residuals
% res   = [di(:,2),-di(:,1)].*(x0-ones(N,1)*x');
Omega = 0;
for n=1:N
    Omega=Omega+(x0(n,:)-x')*w(n)*(eye(2)-di(n,:)'*di(n,:))*(x0(n,:)-x')';
end

if R > 1
    sigma_0 = sqrt(Omega/(N-2));
end
xe=sugr_Point_2D(x,Cxx);
