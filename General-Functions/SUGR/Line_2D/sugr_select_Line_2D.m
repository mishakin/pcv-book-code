%% Select n-th line form list of lines
%
% l = sugr_select_Line_2D(ls,n)
%
% * ls   2D lines, array of structs {l.h(N,:),l.Crr(N,:,:)}
% * n    index of line to be selected
% 
% * l    selected 2D line, struct {l.h,l.Crr}
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_select_Point_2D

function l = sugr_select_Line_2D(ls,n)

l.h    = ls.h(n,:)';
l.Crr  = squeeze(ls.Crr(n,:,:));
l.Jr   = null(l.h');
l.type = 2;

