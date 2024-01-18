%% sugr_ConicM_ConicV transforms conic matrix to conic vector
%
% c = sugr_ConicM_ConicV(C)
% 
% C = 3x3 symmetric matrix
% c = 6-vector
%
% wf 11/2012

function c = sugr_ConicM_ConicV(C)

c = [C(1,1) C(1,2) C(2,2) C(1,3) C(2,3) C(3,3)]';

end