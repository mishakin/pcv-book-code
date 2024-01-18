%% Determine CbRbR from Cee for given Relative Orientaion
%
% [EbR,CbRbR,JebR] = sugr_get_CbRbR_JebR_from_b_R_Cee(b,R,Cee)
%
% input
%    b      = 3x1 basis
%    R      = 3 x 3 rotation matrix
%    Cee    = 9 x 9 Cov(e)
%
% output
%    EbR    = 3 x 3 E-matrix
%    CbRbR  = 5 x 5 CovM of [b_r;r] 
%    JebR   = 9 x 5 Jacobian de/dbR
%
% Wolfgang Förstner 11/2017
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix
function [EbR,CbRbR,JebR] = sugr_get_CbRbR_JebR_from_b_R_Cee(b,R,Cee)

% essential matrix
EbR = calc_S(b)*R';

% Jacobian de/d[b_r;r]
% E  = S(b)*R' ... % rows ri and ei of R and E
% de = vec(dS(b)*R' + S(b)*dR')  (10.339)
% de = [S'(r1); S'(r2); S'(r3)] db + vec(S(b)*R'*S'(dr))
% de = -[S(r1); S(r2); S(r3)] J_r(b) db_r + vec[e1'*S'(dr);e2'*S'(dr);e3'*S'(dr)])
% de = -[S(r1); S(r2); S(r3)] J_r(b) db_r + vec[dr'*S(e1);dr'*S(e3);dr'*S(e3)] 
% de = -[S(r1); S(r2); S(r3)] J_r(b) db_r + vec[dr'*[S(e1);S(e2);S(e3)]]
% using T_mn for permuting the rows of the Jacobian, s. Fackler (2005)
% Notes on matrix calculus -- see Maple jacobian_vecESr_wrt_r.mws.
% de = -[S(r1); S(r2); S(r3)] J_r(b) db_r + T_mn*vec[S'(e1),S'(e2),S'(e3)]*dr
%
% de = JebR * dbR -> dbR = JebR^+ de

% Jacobian wrt b_r: 
Jb = -[calc_S(R(1,:));...
       calc_S(R(2,:));...
       calc_S(R(3,:))]...
     * null(b');          % 9 x 2

 % Jacobian wrt r
z3=zeros(3,1);
Jp = -[ z3,      -EbR(:,3)  EbR(:,2);...
        EbR(:,3) z3        -EbR(:,1);...
       -EbR(:,2)  EbR(:,1)  z3];        % 9 x 3
   
% joint Jacobian wrt [b_r;r]
JebR = [Jb,Jp];                  % 9 x 5

%%% Check numerically
%[~,~,J]=var_prop_classical(@e_from_bR,zeros(5,1),eye(5),[b,R]);

JebRp = (JebR'*JebR)\JebR';      % 5 x 9

% CovM of e=vecE
CbRbR  = JebRp * Cee * JebRp';   % 5 x 5


