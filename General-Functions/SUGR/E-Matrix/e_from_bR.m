%% get e (= vec E) from [db_r,dr]
%
% function for checking Jacobian 9x5 of d(vec E)/d([b_r;r])
%       potentially used in sugr_get_Cee_JebR_from_b_R_CbRbR.m, l. 46
%
% Usage:
%    e = e_from_bR(dbR,bR)
%        db_r = reduced coordinates of b (differential)
%        dr   = differential correction of rotation matrix 
%
%    dbR = 5x1 vector od differentials[db_r;dr] 
%    bR  = 3x4 matrix [b,R], approximate values 
% 
%    e   = 9x1 vector of vec(S(b+db)*(Rot(dr)*R)
%
% Wolfgang Förstner 
% wfoerstn@uni-bonn.de 

function  e = e_from_bR(dbR,bR)

% E  = S(b)*R'
% take b,R as approximate values
b = bR(:,1);
R = bR(:,2:4);
% provide Jacobian Jr
Jrb = null(b');
% determine E
E = calc_S(b+Jrb*dbR(1:2))*(calc_Rot_r(dbR(3:5)/2)*R)';
% finally vectorize
e = E(:);




