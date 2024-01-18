%% Generate N point pairs for given Projection P
%
% [X,y] = sugr_generate_true_2D_point_pairs_Projection_3D_2D(P,N,dX,dY,dZ,boolean_r);
%
% P         = 3x4 matrix
% N         = number of points
% dX,dY,dZ  = size of box around [0,0,0]
% boolean_r = boolean: points should sit random
%                else: points sit in a square (N should be square)
%
% [X,y] = point pairs
%      X.e,y.e      = N x 3, N x 2 matrices of Cartesian point coordinates
%      X.Cee,y.Cee  = N x 3 x 3, N x 2 x 2 with CovM of point pairs
%
% Wolfgang Förstner 02/2013
% wfoerstn@uni-bonn.de
%
% See also sugr_estimation_algebraic_Projection_3D_2D_from_point_pairs
% sugr_estimation_ml_Projection_3D_2D_from_point_pairs

function [X,y] = sugr_generate_true_2D_point_pairs_Projection_3D_2D(P,N,dX,dY,dZ,br)

if br
    for n=1:N
        Xe             = [rand(1)*dX;...
            rand(1)*dY;...
            rand(1)*dZ;...
            1];
        ye             = P * Xe;
        X.e(n,:)       = Xe(1:3)'/Xe(4);
        X.Cee(n,:,:)   = zeros(3);
        y.e(n,:)       = ye(1:2)'/ye(3);
        y.Cee(n,:,:)   = zeros(2);
    end
else
    M = ceil(N^(1/3));
    n=0;
    for k=1:M
        for l=1:M
            for m=1:M
                n = n+1;
                if n <= N
                    Xe = [(k-(M+1)/2)/(M-1);...
                        (l-(M+1)/2)/(M-1);...
                        (m-(M+1)/2)/(M-1);...
                        1];
                    ye           = P * Xe;
                    X.e(n,:)       = Xe(1:3)'/Xe(4);
                    X.Cee(n,:,:)   = zeros(3);
                    y.e(n,:)       = ye(1:2)'/ye(3);
                    y.Cee(n,:,:)   = zeros(2);
                end
            end
        end
    end
end
