%% sugr_ConicV_ConicM transforms conic vector to conic matrix 
%
% c = sugr_ConicV_ConicM(C)
% 
% c = 6-vector
% C = 3x3 symmetric matrix
%
% wf 12/2012

function C = sugr_ConicV_ConicM(c)

C = [c(1) c(2) c(4) ;...
     c(2) c(3) c(5) ;...
     c(4) c(5) c(6)];

end