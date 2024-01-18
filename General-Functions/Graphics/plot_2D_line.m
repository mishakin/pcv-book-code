%% plots a 2D line from homogeneous vector h [3x1]
% 
% Usage:
%   plot_2D_line(h) 
%        plots line h into current figure, using the max extend of
%        current axis object
%   plot_2D_line(h,x1,x2) 
%        plots line h between points x1 and x2, which are
%        either eucledian [2x1] or homogeneous [3x1]
%   h = plot_2D_line(h,varargin) 
%        returns line handle 
%
% 06/17 Susanne Wenzel
% wenzel@igg.uni-bonn.de

function varargout = plot_2D_line(h,varargin)

N = @(x)x./x(3);

h = h./norm(h(1:2));

% get the current axis objects extend
x_lim = get(gca,'xlim');
y_lim = get(gca,'ylim');
       
if nargin<3
    
    % generate standardline, footpoint or origin [+d, -d]
    d = norm([max(abs(x_lim)); max(abs(y_lim))]);

    S3 = [0 -1 0; 1 0 0; 0 0 0];
    % footpoint
    zl0 = N(calc_S(h)*S3*h);
    % direction
    phi = atan2(h(2),h(1));
    % start- and endpoint
    xs = zl0 + d*[cos(phi+pi/2);sin(phi+pi/2);0];
    xe = zl0 + d*[cos(phi-pi/2);sin(phi-pi/2);0];

else
    % given startpoint
    xs = varargin{1};
    if numel(xs)<3, xs = N([xs;1]);end
    % footpoint onto the line, as startpoint for ploting
    xs = calc_S(h)* calc_S(xs)*diag([1,1,0])*h;
    
    % given endpoint
    xe = varargin{2};
    if numel(xe)<3, xe = N([xe;1]);end
    % footpoint onto the line, as endpoint for ploting
    xe = calc_S(h)* calc_S(xe)*diag([1,1,0])*h;

end

handle = plot([xs(1),xe(1)],[xs(2),xe(2)]);
set(gca,'xlim',x_lim);
set(gca,'ylim',y_lim);

if nargout>0
    varargout{1} = handle;
end

