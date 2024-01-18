%% Select n-th point form list of points
%
% x = sugr_select_Point_2D(xs,n)
%
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de
% last change: 9/16 sw
%
% See also sugr_Point_2D


function x = sugr_select_Point_2D(xs,n)

x.h    = xs.h(n,:)';
x.Crr  = squeeze(xs.Crr(n,:,:));
x.type = 1;
