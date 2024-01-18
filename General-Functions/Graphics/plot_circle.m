%% plot circle.
%
% Usage:
%    plot_circle( xo,yo,r, fmt)
%
% xo,yo: centre
%     r: radius
%   fmt: format string
%
% author J. meidow, FGAN-FOM
%
% 2004, Oct:  initial code
% 2005, Jan.: varargs and doc. added
function h = plot_circle(xo,yo,r, fmt, varargin)

if r > 0.0
    h = plot_ellipse([xo,yo,r,r,0.0],fmt);
    if nargin>4
        set(h,varargin{:});
    end
end
