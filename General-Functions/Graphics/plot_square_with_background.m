%% plot box with background
%
% Usage:
%   plot_square_with_background(x,y,w);
%
% x,y:  position
% w: sidelength

function plot_square_with_background(x,y,w,f1,f2)

if nargin < 4
    f1 = 4;
    f2 = 2;
end

hi = w/2;
ho = hi;

plot([x-hi,x-hi,x+hi,x+hi,x-hi],[y-hi,y+hi,y+hi,y-hi,y-hi],'-k','Linewidth',f1);
plot([x-ho,x-ho,x+ho,x+ho,x-ho],[y-ho,y+ho,y+ho,y-ho,y-ho],'-y','Linewidth',f2);
