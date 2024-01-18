%% plot_ellipse: plot 2D ellipse
%
% usage:
% h = plot_ellipse(p,fmt) using parameters p = [xo,yo,a,b,phi]
% h = plot_ellipse(C,fmt) using 3x3 conic C with det(Chh)!>0, Chh = C(1:2,1:2)
% h = plot_ellipse(Cov,fmt) using Cov.mean and Cov.cov covariance matrix
%
% arguments:
%  p - double 5x1, parameters of ellipse p = [xo,yo,a,b,phi]
%                  xo,yo:  centre
%                  a, b:  semi axes
%                  phi:  direction
%  C - double 3x3, conic
%  Cov - struct with fields
%        .mean - double 2x1
%        .cov  - double 2x2 covariance matrix to plot at mean
%  fmt:  optional format string
%
% output:
%      h:  handle to line object
%
% author: J. Meidow
%
% 2004, Oct.: adaptive number of supporting points 
% 2005, Mar.: check input args
% 2016, Sep.: added conic and covariance matrix as input 
%             Susanne Wenzel, wenzel@igg.uni-bonn.de
function h = plot_ellipse(ellipse,fmt)

if nargin<1
    error('not enough input arguments')
end

if nargin>1 && ~ischar(fmt)
    error('2nd input must be string')
end

if nargin==1
    fmt = 'k-';
end

if isstruct(ellipse)
    p(1) = ellipse.mean(1);       % x0
    p(2) = ellipse.mean(2);       % y0
    [~,D,R] = svd(ellipse.cov);
    p(3) = sqrt(D(1,1));        % a
    p(4) = sqrt(D(2,2));        % b
    p(5) = atan2(R(2,1),R(1,1));  % phi
elseif ~isvector(ellipse)
    p = ellipseConic2param(ellipse);
else
    p = ellipse;
end
  
xo = p(1);
yo = p(2);
a = p(3);
b = p(4);
phi = p(5);
    

if a<=0.0
    warning('plot ellipse: semi-major axis a = %g <=0', a)
end
if b<=0.0
    warning('plot ellipse: semi-minor axis b = %g <=0', b)
end

% ds = 10;              % pixel
% perim = 2*(a+b);        % perimeter
% N = round(perim/ds);  % number of supporting points 16...200
% N = min( N, 200);
% N = max( N, 200);
N = 20000;

t = linspace( 0.0, 2.0*pi, N );

u = a*sin(t);
v = b*cos(t);
x = cos(phi)*u -sin(phi)*v +xo;
y = sin(phi)*u +cos(phi)*v +yo;

h = plot(x,y,fmt);

