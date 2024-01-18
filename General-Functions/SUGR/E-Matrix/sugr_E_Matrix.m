%% create E-matrix
%
% Usage
%    E = sugr_E_Matrix(E)         
%       3x3 matrix, needs not be an essential matrix, CovM = 0
%    E = sugr_E_Matrix(b,R)       
%       3x1 vector b, 3x3 rotation R, CovM = 0
%    E = sugr_E_Matrix(b,R,CbRbR) 
%       3x1 vector b, 3x3 rotation R, 5x5 CovM of d_r = [db_r,dr]
%    E = sugr_E_Matrix(b,R,Cee)   
%       3x1 vector b, 3x3 rotation R, 9x9 CovM of d_e = d_vecE,
%       needs  not have the correct rank
%
% output
% E-Matrix = structure
%    .E     = 3 x 3-Matrix, 
%             Frobenius norm = 2
%             det(U)=det(V), for [U,D,V]=svd(E)% 
%    .bR    = 3 x 4 matrix [b,R], if input b, R  
%    .Cee   = 9x9 covariance matrix of e=vecE, rank=5
%    .CbRbR = 5 x 5 covariance matrix of minimal parameters r
%             r(5) = [db_r; dm] for [b,R], if input b, R
%    .JebR  = 8 x 5 Jacobian de/dbR
%    .type  = 25 -> Essential matrix
%
% Wolfgang Förstner 11/2017
% wfoerstn@uni-bonn.de

function E = sugr_E_Matrix(a1,a2,a3)

% default
Cee   = zeros(9);                                                          %#ok<PREALL>
CbRbR = zeros(5);
%
switch nargin
    case 1 % E-matrix, CovM = 0
        [U,~,V] = svd(a1);
        Em      = U*diag([1,1,0])*V'*det(U)*det(V);
        % collect
        E.E     = Em;
        E.CbRbR = CbRbR;
    case 2 % b, R
        b        = a1/norm(a1);
        E.E      = calc_S(b)*a2';
        E.bR     = [b,a2];
        E.CbRbR  = CbRbR;
    case 3
        r=size(a3,1);
        switch r
            case 5 % basis, rotation, CbRbR              
                factor  = 1/norm(a1);
                b       = a1*factor;
                CbRbR   = a3*factor^2;
                [EbR,Cee,JebR]=...
                    sugr_get_Cee_JebR_from_b_R_CbRbR(b,a2,CbRbR);
                % normalize E
                [U,~,V] = svd(EbR);
                Em      = U*diag([1,1,0])*V'*det(U)*det(V);
                % collect
                E.E     = Em;
                E.bR    = [b,a2];
                E.Cee   = Cee;     
                E.CbRbR = CbRbR;
                E.JebR  = JebR;
            case 9 % basis, rotation, Cee
                factor  = 1/norm(a1);
                b       = a1*factor;
                Cee     = a3*factor^2;
                [EbR,CbRbR,JebR]=...
                    sugr_get_CbRbR_JebR_from_b_R_Cee(b,a2,Cee); 
                % normalize E
                [U,~,V] = svd(EbR);
                Em      = U*diag([1,1,0])*V'*det(U)*det(V);             
                % normalize E and get JebR
                E.E     = Em;
                E.bR    = [b,a2];
                E.CbRbR = CbRbR;
                E.Cee   = JebR*E.CbRbR*JebR';
                E.JebR  = JebR;
        end
end
E.type = 25;
