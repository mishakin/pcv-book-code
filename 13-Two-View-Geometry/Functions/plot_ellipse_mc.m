% plot ellipse from mean and covariance matrix
%
% plot_ellipse_mc(mu,cov,p_typ)
% mu = centre
% cov = covariance matrix
% p_typ = point type
%
% Wolfgang Förstner 10/2010
% wfoerstn@uni-bonn.de


function h = plot_ellipse_mc(mu,cov,p_typ,lw)

if nargin < 4
    lw = 2;
end

e = sqrt(eigs(cov));
a = e(1);
b = e(2);
theta = atan2(2*cov(1,2),cov(1,1)-cov(2,2))/2;

np = 100;
ang = (0:np)*2*pi/np;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = mu*ones(size(ang)) + R*[cos(ang)*a; sin(ang)*b];

h = plot( pts(1,:), pts(2,:) , p_typ,'LineWidth',lw);
