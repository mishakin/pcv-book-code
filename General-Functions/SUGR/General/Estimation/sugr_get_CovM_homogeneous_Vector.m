%% get CovM of homogeneous vector x.h
%
% Chh = sugr_get_CovM_homogeneous_Vector(x); 
% 
% * x =  struct, homogeneous object, minimal representation {x.h,x.Crr}
%
% * Chh = full cov. matrix of homogeneous element x.h

% see PCV (10.28)
%
% Wolfgang Förstner 2016-09-02
% wfoerstn@uni-bonn.de 

function Chh = sugr_get_CovM_homogeneous_Vector(x)

J = null(x.h');        % d x (d-1)

Chh = J * x.Crr * J';

