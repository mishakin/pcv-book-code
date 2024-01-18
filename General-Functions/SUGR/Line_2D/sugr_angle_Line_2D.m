%% Get angle between two 2D lines 
%
% angle = sugr_angle_Line_2D(l,m)
%
% * l, m   2D line, struct {l.h,l.Crr}
%
% * angle.a  = angle between lines, from l to m in [0,2*pi)%
%        .sa = standard deviation angle negativ if one line is infinte%
%        .dof = degrees of freedom = 1
%
% Wolfgang Förstner 2/2011
% wfoerstn@uni-bonn.de 
%
% See also sugr_Line_2D 

function angle = sugr_angle_Line_2D(l,m)


if sugr_get_isfinite_Line_2D(l) && sugr_get_isfinite_Line_2D(m)

    % Hessean parameters of l and m
    [le, Clee] = sugr_get_Euclidean_Line_2D(l);
    [me, Cmee] = sugr_get_Euclidean_Line_2D(m);

    % angle and standard deviation
    angle.a = mod(me(1) - le(1),2*pi);
    angle.sa = sqrt(Clee(1,1) + Cmee(1,1));

else        % one line at infinity
    angle.a = 0;
    angle.sa = -1;
end
angle.dof = 1;







